import scipy.io as scio
import scipy.sparse as scsp
import scipy.sparse.linalg as scspli
from dolfin import *
from linearized_aux import u_uncompress
from linearized_aux import u_compress
import numpy as np


# set floating point error handling
np.seterr(all="raise", divide="raise", over="raise", under="raise", invalid="raise")


class LinearizedCtrl():
    def __init__(self, ref, RE, pertubationeps, dt, T):

        # parameters
        self.ref = ref
        self.RE = RE
        self.pertubationeps = pertubationeps
        self.dt = dt
        self.T = T
        self.t = 0
        self.k = 0

        self.logv = np.zeros((np.ceil(self.T/self.dt)+1, 6))

        # mesh and function spaces
        self.mesh = Mesh(const.MESH_XML(ref))
        self.V = VectorFunctionSpace(self.mesh, const.LINEARIZED_SIM_V, const.LINEARIZED_SIM_V_DIM)
        self.Q = FunctionSpace(self.mesh, const.LINEARIZED_SIM_Q, const.LINEARIZED_SIM_Q_DIM)
        self.u_stat = Function(self.V, const.STATIONARY_U_XML(ref, RE))

        # compress system
        self.Mcps = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_M_MTX(self.ref, self.RE))
        self.Mlowercps = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_MLOWER_MTX(self.ref, self.RE))
        self.Muppercps = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_MUPPER_MTX(self.ref, self.RE))
        self.Scps = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_S_MTX(self.ref, self.RE))
        self.Rcps = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_R_MTX(self.ref, self.RE))
        self.Kcps = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_K_MTX(self.ref, self.RE))
        self.Gcps = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_G_MTX(self.ref, self.RE))
        self.GTcps = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_GT_MTX(self.ref, self.RE))
        self.Bcps = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_B_MTX(self.ref, self.RE))
        self.Ccps = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_C_MTX(self.ref, self.RE))
        self.inner_nodes = np.loadtxt(const.ASSEMBLER_COMPRESS_CTRL_INNER_NODES(self.ref, self.RE), dtype=np.int64)

        # system sizes
        self.nv, self.np = self.Gcps.shape

        # set other attributes
        # store feedback matrix
        self.Kinfcps = scio.mmread(const.LINEARIZED_CTRL_KINF_CPS_MTX(self.ref, self.RE))

        # lift input matrix and feedback
        self.Kinfsys = np.vstack((self.Kinfcps, np.zeros((self.np, self.Kinfcps.shape[1]))))
        self.Bsys = -1*self.dt * np.vstack((self.Bcps, np.zeros((self.np, self.Bcps.shape[1]))))

        # visualisation
        self.u_dolfin = Function(self.V)
        self.u_file = File(const.LINEARIZED_CTRL_U_PVD(self.ref, self.RE))
        self.udelta_dolfin = Function(self.V)
        self.udelta_file = File(const.LINEARIZED_CTRL_U_DELTA_PVD(self.ref, self.RE))

        # define state
        self.uk_sys = self.u_stat.vector().array().reshape(len(self.u_stat.vector().array()), )[self.inner_nodes]
        self.uk_sys *= self.pertubationeps

        # build system matrices
        self.Asys = -self.Scps-self.Rcps-self.Kcps-self.Muppercps-self.Mlowercps
        u = scsp.hstack([self.Mcps - self.dt*self.Asys, self.dt*(-self.Gcps)])
        l = scsp.hstack([self.dt * (-self.GTcps), scsp.csr_matrix((self.np, self.np))])
        self.Msys_ode = scsp.vstack([u, l]).tocsc()
        self.Msys_lift = scsp.vstack([self.Mcps, scsp.csr_matrix((self.np, self.nv))])

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

        if self.t < const.LINEARIZED_CTRL_START_CONTROLLING:
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
        self.logv[self.k, 0] = self.t
        self.logv[self.k, 1] = np.linalg.norm(self.uk_sys)

        # log control
        uc = self.Kinfcps.T*self.uk_sys
        self.logv[self.k, 2] = uc[0, 0]
        self.logv[self.k, 3] = uc[0, 1]

        # log output
        uout = self.Ccps*self.uk_sys
        self.logv[self.k, 4] = uout[0, 0]
        self.logv[self.k, 5] = uout[0, 1]


    def _save(self):

        if (self.k-1) % const.LINEARIZED_CTRL_SAVE_FREQ == 0:
            # uncompress and save pvd for udelta
            uk_uncps = u_uncompress(self.V, self.uk_sys, self.inner_nodes)
            self.udelta_dolfin.vector().set_local(uk_uncps)
            self.udelta_file << self.udelta_dolfin

            # add stationary and save pvd
            v = uk_uncps + self.u_stat.vector().array()
            self.u_dolfin.vector().set_local(v)
            self.u_file << self.u_dolfin

    def solve_ode(self):
        while self.t < (self.T + DOLFIN_EPS):
            self.log()

            # print info
            if self.k % int(const.LINEARIZED_CTRL_INFO*(self.T/self.dt)) == 0:
                print "{0:.2f}%\t t={1:.3f}\t ||u_delta||={2:e}".format(self.t/self.T*100, self.logv[self.k, 0], self.logv[self.k, 1])

            self._save()
            self.uk_next()

    def save_log(self):
        np.savetxt(const.LINEARIZED_CTRL_LOG(self.ref, self.RE), self.logv)