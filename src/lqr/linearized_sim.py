import scipy.io as scio
import scipy.sparse as scsp
import scipy.sparse.linalg as scspli
import numpy as np
from dolfin import *
from src.aux import profile, gettime
#set floating point error handling
np.seterr(all="raise", divide="raise", over="raise", under="raise", invalid="raise")


class LinearizedSim():
    @profile
    def __init__(self, const, ref, RE):

        # parameters
        self.ref = ref
        self.RE = RE
        self.const = const
        self.pertubationeps = self.const.LINEARIZED_SIM_PERTUBATIONEPS
        self.dt = self.const.LINEARIZED_SIM_DT
        self.T = self.const.LINEARIZED_SIM_T
        self.t = 0.0
        self.k = 0
        self.save_freq = int(1.0/(self.dt*self.const.LINEARIZED_SIM_SAVE_PER_S))

        # log
        self.logv = np.zeros((int(self.T/(self.dt*self.save_freq)+1), 2))
        self.klog = 0

        # mesh and function spaces
        self.mesh = Mesh(const.MESH_XML(ref))
        self.V = VectorFunctionSpace(self.mesh, const.V, const.V_DIM)
        self.Q = FunctionSpace(self.mesh, const.Q, const.Q_DIM)
        self.u_stat = Function(self.V, const.STATIONARY_U_XML(ref, RE))

        # read compress system for simuation
        names = ["M", "M_BOUNDARY_CTRL", "S", "R", "K", "G", "GT", "B", "C"]
        self.mat = {}
        for name in names:
            self.mat[name] = scio.mmread(const.ASSEMBLER_COMPRESS_SIM_NAME_MTX(ref, name, RE))
        self.inner_nodes = np.loadtxt(const.ASSEMBLER_COMPRESS_SIM_INNERNODES_DAT(ref, RE), dtype=np.int64)

        # system sizes
        self.nv, self.np = self.mat["G"].shape
        self.ninner = self.inner_nodes.size

        # visualization
        self.u_dolfin = Function(self.V)
        self.u_dolfin.rename(*self.const.PVD_U_LABEL_NAME)
        self.u_file = File(const.LINEARIZED_SIM_U_PVD(self.ref, self.RE))
        self.udelta_dolfin = Function(self.V)
        self.udelta_dolfin.rename(*self.const.PVD_U_LABEL_NAME)
        self.udelta_file = File(const.LINEARIZED_SIM_U_DELTA_PVD(self.ref, self.RE))

        # define state and uncompress state
        self.uk_sys = np.take(self.u_stat.vector().array(), self.inner_nodes)
        self.uk_sys *= self.pertubationeps
        self.uk_uncps = np.zeros((self.V.dim(),))

        # build system matrices
        self.Asys = -self.mat["S"]-self.mat["R"]-self.mat["K"]
        u = scsp.hstack([self.mat["M"] - self.dt*self.Asys, self.dt*(-self.mat["G"])])
        l = scsp.hstack([self.dt * (-self.mat["GT"]), scsp.csr_matrix((self.np, self.np))])
        self.Msys_ode = scsp.vstack([u, l]).tocsc()
        self.Msys_lift = scsp.vstack([self.mat["M"], scsp.csr_matrix((self.np, self.nv))]).tocsr()

        # build solver
        # warning for higher refinements there are sometimes error change these parameters or ilu<->lu
        self.Msys_solver = scspli.splu(self.Msys_ode).solve

        # test function and form for assembleN
        self.w_test = TestFunction(self.V)
        self.N = inner(self.w_test, grad(self.u_dolfin) * self.u_dolfin)*dx
        self.N_uk_uk = np.zeros((self.np+self.ninner,))

    @profile
    def assembleN(self,uk):

        # insert values at inner nodes
        self.uk_uncps[self.inner_nodes] = uk

        # update vector, assemble form  and compress vector
        self.u_dolfin.vector().set_local(self.uk_uncps)
        Nassemble = assemble(self.N)

        # compress to inner nodes and fill into upper block
        self.N_uk_uk[:self.ninner] = np.take(Nassemble.array(), self.inner_nodes)

    @profile
    def uk_next(self):
        # compute new rhs (prediction)
        self.assembleN(self.uk_sys)
        Msys_lift_uk = self.Msys_lift * self.uk_sys
        rhs = Msys_lift_uk - self.dt*self.N_uk_uk
        self.uk_sys = self.Msys_solver(rhs)[0:self.nv]
        # try to correct solution
        self.correction(Msys_lift_uk)
        self.k += 1
        self.t += self.dt

    def correction(self, Msys_lift_uk):
        # correct predicted solution

        for i in range(self.const.LINEARIZED_SIM_CORRECTION_STEPS):
            # assemble nonlinear right hand side term
            self.assembleN(self.uk_sys)
            ukpk_sys_correction = self.Msys_solver(-self.dt*self.N_uk_uk+Msys_lift_uk)
            uk_sys_correction = ukpk_sys_correction[0:self.nv]

            # compute residual for implicit euler scheme
            self.assembleN(uk_sys_correction)
            res = np.linalg.norm(self.Msys_ode*ukpk_sys_correction+self.dt*self.N_uk_uk-Msys_lift_uk)

            # show info
            if self.k % int(self.const.LINEARIZED_SIM_INFO*(self.T/self.dt)) == 0:
                diff = np.linalg.norm(self.uk_sys-uk_sys_correction)
                print "Correction Step {0:d}\t||uk_sys-uk_sys_correction||={1:e}\t||res (implicit euler)||={2:e}".format(i, diff, res)

            self.uk_sys = uk_sys_correction
            if res < self.const.LINEARIZED_SIM_CORRECTION_RES:
                break

    def log(self):
        if self.k % self.save_freq == 0:
            self.logv[self.klog, 0] = self.t
            self.logv[self.klog, 1] = np.linalg.norm(self.uk_sys)
            self.klog += 1

    def save(self):
        if self.k % self.save_freq == 0:

            # uncompress and save pvd for udelta
            self.uk_uncps[self.inner_nodes] = self.uk_sys
            self.udelta_dolfin.vector().set_local(self.uk_uncps)
            self.udelta_file << (self.udelta_dolfin, self.t)

            # add stationary and save pvd
            v = self.uk_uncps + self.u_stat.vector().array()
            self.u_dolfin.vector().set_local(v)
            self.u_file << (self.u_dolfin, self.t)

    @profile
    def solve_ode(self):
        while self.t < (self.T + DOLFIN_EPS):
            self.save()
            self.log()

            # print info
            if self.k % int(self.const.LINEARIZED_SIM_INFO*(self.T/self.dt)) == 0:
                print gettime(), "{0:.2f}%\t t={1:.3f}\t ||u_delta||={2:e}".format(self.t/self.T*100, self.t, np.linalg.norm(self.uk_sys))

            self.uk_next()

    def save_log(self):
        # self.logv = self.logv[:self.k, :]
        np.savetxt(self.const.LINEARIZED_SIM_LOG(self.ref, self.RE), self.logv)

