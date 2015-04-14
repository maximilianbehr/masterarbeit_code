import scipy.io as scio
import scipy.sparse as scsp
import scipy.sparse.linalg as scspli
from dolfin import *
import numpy as np
from src.aux import profile, gettime
from linearized_aux import stable_timestep


#set floating point error handling
np.seterr(all="raise", divide="raise", over="raise", under="raise", invalid="raise")

class LinearizedSim():

    def __init__(self, const, ref, RE):

        # parameters
        self.ref = ref
        self.RE = RE
        self.const = const
        self.pertubationeps = self.const.LINEARIZED_SIM_PERTUBATIONEPS
        self.dt = self.const.LINEARIZED_SIM_DT
        self.T = self.const.LINEARIZED_SIM_T

        # mesh and function spaces
        self.mesh = Mesh(const.MESH_XML(ref))
        self.V = VectorFunctionSpace(self.mesh, const.V, const.V_DIM)
        self.Q = FunctionSpace(self.mesh, const.Q, const.Q_DIM)
        self.u_stat = Function(self.V, const.STATIONARY_U_XML(ref, RE))

        if const.LINEARIZED_SIM_STABLE_DT:
            self.dt = stable_timestep(self.T, self.const.GET_NU_FLOAT(self.RE), self.const.U, self.mesh.hmin())

        print "dt = {0:e}".format(self.dt)

        self.t = 0.0
        self.k = 0
        self.save_freq = int(1.0/(self.dt*self.const.LINEARIZED_SIM_SAVE_PER_S))

        # log
        self.logv = np.zeros((int(self.T/(self.dt*self.save_freq)+1), 2))
        self.klog = 0
        self.kinfo = int(self.const.LINEARIZED_SIM_INFO*(self.T/self.dt))

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
        u_pertubated = project(self.pertubationeps*self.u_stat, self.V)
        self.uk_sys = np.take(u_pertubated.vector().array(), self.inner_nodes)
        self.uk_uncps = np.zeros((self.V.dim(),))
        self.uk_uncps[self.inner_nodes] = self.uk_sys

        # build system matrices
        self.Asys = -self.mat["S"]-self.mat["R"]-self.mat["K"]
        u = scsp.hstack([self.mat["M"] - self.dt*self.Asys, -self.dt*(self.mat["G"])])
        l = scsp.hstack([-self.dt * (self.mat["GT"]), scsp.csr_matrix((self.np, self.np))])
        self.Msys_ode = scsp.vstack([u, l]).tocsc()
        self.Msys_lift = scsp.vstack([self.mat["M"], scsp.csr_matrix((self.np, self.nv))]).tocsr()

        # build solver
        self.Msys_solver = scspli.splu(self.Msys_ode.tocsc()).solve

        # test function and form for assembleN
        self.w_test = TestFunction(self.V)
        self.N = inner(self.w_test, grad(self.u_dolfin) * self.u_dolfin)*dx
        self.N_uk_uk = np.zeros((self.np+self.ninner,))
        self.assembleN(self.uk_sys)

    def assembleN(self, uk):

        # insert values at inner nodes
        self.uk_uncps[self.inner_nodes] = uk

        # update vector, assemble form  and compress vector
        self.u_dolfin.vector().set_local(self.uk_uncps)
        Nassemble = assemble(self.N)

        # compress to inner nodes and fill into upper block
        self.N_uk_uk[:self.ninner] = np.take(Nassemble.array(), self.inner_nodes, axis=0)


    def uk_next(self):
        # compute new rhs (prediction)
        self.assembleN(self.uk_sys)
        Msys_lift_uk = self.Msys_lift * self.uk_sys
        self.uk_sys = self.Msys_solver(Msys_lift_uk - self.dt*self.N_uk_uk)[0:self.ninner]
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
            self.uk_sys = uk_sys_correction

            # compute residual for implicit euler scheme only in every third step to improve perfomance
            if (i % self.const.LINEARIZED_SIM_CORRECTION_RES_MOD) == 0:
                self.assembleN(uk_sys_correction)
                res = np.linalg.norm(self.Msys_ode*ukpk_sys_correction+self.dt*self.N_uk_uk-Msys_lift_uk)

                # show info
                if self.k % self.kinfo == 0:
                    print "Correction Step {0:d}\t||res (implicit euler)||={1:e}".format(i,  res)

                if np.isnan(res):
                    raise ValueError('nan during computation')

                if res < self.const.LINEARIZED_SIM_CORRECTION_RES:
                    break


    def log(self):
        if self.k % self.save_freq == 0:
            # nrm = np.linalg.norm(self.uk_sys)
            self.uk_uncps[self.inner_nodes] = self.uk_sys
            self.udelta_dolfin.vector().set_local(self.uk_uncps)
            nrm = norm(self.udelta_dolfin, "L2", mesh=self.mesh)
            self.logv[self.klog, 0] = self.t
            self.logv[self.klog, 1] = nrm

            self.klog += 1

            if np.isnan(nrm):
                print "nan during computation"
                return False

            if nrm < self.const.LINEARIZED_SIM_RES:
                print "convergenced reached"
                return False

        return True


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

    def solve_ode(self):
        while self.t < (self.T + DOLFIN_EPS):
            self.save()
            if not self.log():
                break

            # print info
            if self.k % self.kinfo == 0:
                self.uk_uncps[self.inner_nodes] = self.uk_sys
                self.udelta_dolfin.vector().set_local(self.uk_uncps)
                nrm = norm(self.udelta_dolfin, "L2", mesh=self.mesh)
                print gettime(), "{0:.2f}%\t t={1:.3f}\t ||u_delta||={2:e}".format(self.t/self.T*100, self.t, nrm)

            self.uk_next()

    def save_log(self):
        # self.logv = self.logv[:self.k, :]
        np.savetxt(self.const.LINEARIZED_SIM_LOG(self.ref, self.RE), self.logv)

