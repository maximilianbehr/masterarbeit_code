import scipy.io as scio
import scipy.sparse as scsp
import scipy.sparse.linalg as scspli
from dolfin import *
import numpy as np


# set floating point error handling
np.seterr(all="raise", divide="raise", over="raise", under="raise", invalid="raise")


class LinearizedCtrl():
    def __init__(self, const, ref, RE):

        # parameters
        self.ref = ref
        self.RE = RE
        self.const = const
        self.dt = self.const.LINEARIZED_CTRL_DT
        self.T = self.const.LINEARIZED_CTRL_T
        self.pertubationeps = self.const.LINEARIZED_CTRL_PERTUBATIONEPS
        self.t = 0.0
        self.k = 0
        self.save_freq = int(1.0/(self.dt*self.const.LINEARIZED_SIM_SAVE_PER_S))


        # mesh and function spaces
        self.mesh = Mesh(const.MESH_XML(ref))
        self.V = VectorFunctionSpace(self.mesh, const.V, const.V_DIM)
        self.Q = FunctionSpace(self.mesh, const.Q, const.Q_DIM)
        self.u_stat = Function(self.V, const.STATIONARY_U_XML(ref, RE))

        # read compress system
        names = ["M", "M_BOUNDARY_CTRL", "S", "R", "K", "G", "GT", "B", "C", "Kinf"]
        self.mat = {}
        for name in names:
            self.mat[name] = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_NAME_MTX(ref, name, RE))
        self.inner_nodes = np.loadtxt(const.ASSEMBLER_COMPRESS_CTRL_INNERNODES_DAT(ref, RE), dtype=np.int64)

        # log time, norm(u_delta), outputs, inputs
        self.logv = np.zeros((int(self.T/(self.dt*self.save_freq)+1),\
                              2+self.mat["C"].shape[0]+self.mat["B"].shape[1]))
        self.klog = 0

        # system sizes
        self.nv, self.np = self.mat["G"].shape
        self.ninner = self.inner_nodes.size

        # lift input matrix and feedback
        self.Kinfsys = np.vstack((self.mat["Kinf"], np.zeros((self.np, self.mat["Kinf"].shape[1]))))
        self.Bsys = -1*self.dt * np.vstack((self.mat["B"], np.zeros((self.np, self.mat["B"].shape[1]))))

        # visualisation
        self.u_dolfin = Function(self.V)
        self.u_file = File(const.LINEARIZED_CTRL_U_PVD(self.ref, self.RE))
        self.udelta_dolfin = Function(self.V)
        self.udelta_file = File(const.LINEARIZED_CTRL_U_DELTA_PVD(self.ref, self.RE))

        # define state
        self.uk_sys = np.take(self.u_stat.vector().array(), self.inner_nodes)
        self.uk_sys *= self.pertubationeps
        self.uk_uncps = np.zeros((self.V.dim(),))

        # build system matrices
        self.Asys = -self.mat["S"]-self.mat["R"]-self.mat["K"]-self.mat["M_BOUNDARY_CTRL"]
        u = scsp.hstack([self.mat["M"] - self.dt*self.Asys, self.dt*(-self.mat["G"])])
        l = scsp.hstack([self.dt * (-self.mat["GT"]), scsp.csr_matrix((self.np, self.np))])
        self.Msys_ode = scsp.vstack([u, l]).tocsc()
        self.Msys_lift = scsp.vstack([self.mat["M"], scsp.csr_matrix((self.np, self.nv))]).tocsr()

        # warning for higher refinements there are sometimes error change these parameters or ilu<->lu
        self.Msys_solver = scspli.splu(self.Msys_ode)

        # test function and form for assembleN
        self.w_test = TestFunction(self.V)
        self.N = inner(self.w_test, grad(self.u_dolfin) * self.u_dolfin)*dx
        self.N_uk_uk = np.zeros((self.np+self.ninner,))

    def assembleN(self):

        # insert values at inner nodes
        self.uk_uncps[self.inner_nodes] = self.uk_sys

        # update vector, assemble form  and compress vector
        self.u_dolfin.vector().set_local(self.uk_uncps)
        Nassemble = assemble(self.N)

        # compress to inner nodes and fill into upper block
        self.N_uk_uk[:self.ninner] = np.take(Nassemble.array(), self.inner_nodes)

    def uk_next(self):

        # compute new rhs
        self.assembleN()
        rhs = self.Msys_lift * self.uk_sys - self.dt*self.N_uk_uk

        if self.t < self.const.LINEARIZED_CTRL_START_CONTROLLING:
            # print "no feed"
            temp = self.Msys_solver.solve(rhs)

        else:
            # perform sherman morrison woodbury
            Minvrhs = self.Msys_solver.solve(rhs)
            temp = np.dot(self.Kinfsys.T, Minvrhs)

            Minvb = self.Msys_solver.solve(self.Bsys)
            temp2 = np.dot(self.Kinfsys.T, Minvb)
            temp2 = np.eye(temp2.shape[0]) - temp2
            temp  = np.linalg.solve(temp2, temp)

            temp = np.dot(self.Bsys, temp)
            temp = Minvrhs + self.Msys_solver.solve(temp)

        # lift down vector
        self.uk_sys = temp[0:self.nv]

        self.k += 1
        self.t += self.dt

    def log(self):
        # log time and u delta
        if self.k % self.save_freq == 0:
            self.logv[self.klog, 0] = self.t
            self.logv[self.klog, 1] = np.linalg.norm(self.uk_sys)

            # log control
            uc = np.dot(self.mat["Kinf"].T, self.uk_sys)
            uc.ravel()
            self.logv[self.klog, 2:2+uc.shape[0]] = uc

            # log output
            uout = np.dot(self.mat["C"], self.uk_sys)
            uout.ravel()
            self.logv[self.klog, 4:4+uout.shape[0]] = uout
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

    def solve_ode(self):
        while self.t < (self.T + DOLFIN_EPS):
            self.save()
            self.log()

            # print info
            if self.k % int(self.const.LINEARIZED_CTRL_INFO*(self.T/self.dt)) == 0:
                print "{0:.2f}%\t t={1:.3f}\t ||u_delta||={2:e}".format(self.t/self.T*100, self.t, np.linalg.norm(self.uk_sys))

            self.uk_next()

    def save_log(self):
        #self.logv = self.logv[:self.k, :]
        np.savetxt(self.const.LINEARIZED_CTRL_LOG(self.ref, self.RE), self.logv)