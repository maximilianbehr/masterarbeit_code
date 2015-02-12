from linearized_aux import u_uncompress
from linearized_aux import u_compress
import scipy.io as scio
import scipy.sparse as scsp
import scipy.sparse.linalg as scspli
import numpy as np

from dolfin import *


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
        self.t = 0
        self.k = 0

        # log
        self.logv = np.zeros((np.ceil(self.T/self.dt)+1, 2))

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
        self.nv = self.mat["M"].shape[0]
        self.np = self.mat["G"].shape[1]

        # visualization
        self.u_dolfin = Function(self.V)
        self.u_file = File(const.LINEARIZED_SIM_U_PVD(self.ref, self.RE))
        self.udelta_dolfin = Function(self.V)
        self.udelta_file = File(const.LINEARIZED_SIM_U_DELTA_PVD(self.ref, self.RE))

        # define state
        self.uk_sys = self.u_stat.vector().array().reshape(len(self.u_stat.vector().array()), )[self.inner_nodes]
        self.uk_sys *= self.pertubationeps

        # build system matrices
        self.Asys = -self.mat["S"]-self.mat["R"]-self.mat["K"]
        u = scsp.hstack([self.mat["M"] - self.dt*self.Asys, self.dt*(-self.mat["G"])])
        l = scsp.hstack([self.dt * (-self.mat["GT"]), scsp.csr_matrix((self.np, self.np))])
        self.Msys_ode = scsp.vstack([u, l]).tocsc()
        self.Msys_lift = scsp.vstack([self.mat["M"], scsp.csr_matrix((self.np, self.nv))])

        # build solver
        self.Msys_solver = scspli.spilu(self.Msys_ode)

    def assembleN(self):

        # Function and TestFunction
        w_test = TestFunction(self.V)

        # insert values at inner nodes
        uk = u_uncompress(self.V, self.uk_sys, self.inner_nodes)

        # assemble and compress vector
        self.u_dolfin.vector().set_local(uk)
        N = inner(grad(self.u_dolfin) * self.u_dolfin, w_test)*dx
        Nassemble = assemble(N)
        self.N_uk_uk = u_compress(Nassemble.array().reshape((self.V.dim(),)), self.inner_nodes)

        # lift
        self.N_uk_uk = np.append(self.N_uk_uk, np.zeros((self.np,)))

    def uk_next(self):

        # compute new rhs
        self.assembleN()
        rhs = self.Msys_lift * self.uk_sys - self.dt*self.N_uk_uk

        self.uk_sys = self.Msys_solver.solve(rhs)[0:self.nv]
        self.k += 1
        self.t += self.dt

    def log(self):
        self.logv[self.k, 0] = self.t
        self.logv[self.k, 1] = np.linalg.norm(self.uk_sys)

    def save(self):
        if (self.k-1) % self.const.LINEARIZED_SIM_SAVE_FREQ == 0:

            # uncompress and save pvd for udelta
            uk_uncps = u_uncompress(self.V, self.uk_sys, self.inner_nodes)
            self.udelta_dolfin.vector().set_local(uk_uncps)
            self.udelta_file <<(self.udelta_dolfin, self.t)

            # add stationary and save pvd
            v = uk_uncps + self.u_stat.vector().array()
            self.u_dolfin.vector().set_local(v)
            self.u_file << (self.u_dolfin, self.t)

    def solve_ode(self):
        while self.t < (self.T + DOLFIN_EPS):
            self.log()

            # print info
            if self.k % int(self.const.LINEARIZED_SIM_INFO*(self.T/self.dt)) == 0:
                print "{0:.2f}%\t t={1:.3f}\t ||u_delta||={2:e}".format(self.t/self.T*100, self.logv[self.k, 0], self.logv[self.k, 1])

            self.save()
            self.uk_next()


    def save_log(self):
        self.logv = self.logv[0:self.k, :]
        np.savetxt(self.const.LINEARIZED_SIM_LOG(self.ref, self.RE), self.logv)




