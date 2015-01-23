import src.karman_const as const
from linearized_aux import compress_mat
from linearized_aux import inner_outer_nodes
from linearized_aux import u_uncompress
from linearized_aux import u_compress
from src.mesh.karman import GammaRight
import scipy.io as scio
import scipy.sparse as scsp
import scipy.sparse.linalg as scspli
import numpy as np

from dolfin import *


#set floating point error handling
np.seterr(all="raise", divide="raise", over="raise", under="raise", invalid="raise")


class LinearizedSim():
    def __init__(self, ref, RE, pertubationeps, dt, T):

        # parameters
        self.ref = ref
        self.RE = RE
        self.pertubationeps = pertubationeps
        self.dt = dt
        self.T = T
        self.t = 0
        self.k = 1

        # mesh and function spaces
        self.mesh = Mesh(const.MESH_XML(ref))
        self.V = VectorFunctionSpace(self.mesh, const.LINEARIZED_SIM_V, const.LINEARIZED_SIM_V_DIM)
        self.Q = FunctionSpace(self.mesh, const.LINEARIZED_SIM_Q, const.LINEARIZED_SIM_Q_DIM)
        self.u_stat = Function(self.V, const.STATIONARY_U_XML(ref, RE))

        # define attributes v_inner_nodes and v_outer_nodes
        # GammaRight Boundary Part is added to the inner nodes
        self.inner_nodes, self.outer_nodes = inner_outer_nodes(self.V, [GammaRight()])

        # read matrices
        M = scio.mmread(const.LQR_M_MTX(ref, RE))
        S = scio.mmread(const.LQR_S_MTX(ref, RE))
        R = scio.mmread(const.LQR_R_MTX(ref, RE))
        K = scio.mmread(const.LQR_K_MTX(ref, RE))
        G = scio.mmread(const.LQR_G_MTX(ref, RE))
        GT = scio.mmread(const.LQR_GT_MTX(ref, RE))

        # compress system
        self.Mcps = compress_mat(M, self.inner_nodes)
        self.Scps = compress_mat(S, self.inner_nodes)
        self.Rcps = compress_mat(R, self.inner_nodes)
        self.Kcps = compress_mat(K, self.inner_nodes)
        self.Gcps = compress_mat(G, self.inner_nodes)
        self.GTcps = compress_mat(GT, self.inner_nodes)

        # system sizes
        self.nv = self.Mcps.shape[0]
        self.np = self.Gcps.shape[1]


        # attributes for post init part
        self.u_dolfin = None
        self.udelta_dolfin = None
        self.u_file = None
        self.udelta_file = None
        self.uk_sys = None
        self.Asys = None
        self.Msys_ode = None
        self.Msys_lift = None
        self.Msys_solver = None
        self.N_uk_uk = None

    def postinit(self):

        # visualisation
        self.u_dolfin = Function(self.V)
        self.u_file = File(const.LINEARIZED_SIM_U_PVD(self.ref, self.RE))
        self.udelta_dolfin = Function(self.V)
        self.udelta_file = File(const.LINEARIZED_SIM_U_DELTA_PVD(self.ref, self.RE))

        # define state
        self.uk_sys = self.u_stat.vector().array().reshape(len(self.u_stat.vector().array()), )[self.inner_nodes]
        self.uk_sys *= self.pertubationeps

        # build system matrices
        self.Asys = -self.Scps-self.Rcps-self.Kcps
        u = scsp.hstack([self.Mcps - self.dt*self.Asys, self.dt*(-self.Gcps)])
        l = scsp.hstack([self.dt * (-self.GTcps), scsp.csr_matrix((self.np, self.np))])
        self.Msys_ode = scsp.vstack([u, l]).tocsc()
        self.Msys_lift = scsp.vstack([self.Mcps, scsp.csr_matrix((self.np, self.nv))])

        # build solver
        self.Msys_solver = scspli.spilu(self.Msys_ode)


    def save_compressed_matrices(self):
        """save compress system matrices"""

        matrices = [self.Mcps,
                    self.Scps,
                    self.Rcps,
                    self.Kcps,
                    self.Gcps,
                    self.GTcps]

        files = [const.LINEARIZED_SIM_M_CPS_MTX,
                 const.LINEARIZED_SIM_S_CPS_MTX,
                 const.LINEARIZED_SIM_R_CPS_MTX,
                 const.LINEARIZED_SIM_K_CPS_MTX,
                 const.LINEARIZED_SIM_G_CPS_MTX,
                 const.LINEARIZED_SIM_GT_CPS_MTX]

        for idx in range(len(matrices)):
            with open(files[idx](self.ref, self.RE), "w") as handle:
                scio.mmwrite(handle, matrices[idx])

        with open(const.LINEARIZED_SIM_INNER_NODES(self.ref, self.RE), "w") as handle:
            np.save(handle, self.inner_nodes)

        with open(const.LINEARIZED_SIM_OUTER_NODES(self.ref, self.RE), "w") as handle:
            np.save(handle, self.outer_nodes)


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

        #turn control on
        #Rhs -= (self.Bcompress[:,0]*0.01*sin(self.t)+ self.Bcompress[:,1]*0.02*cos(self.t))

        # solve for the next step
        temp = self.Msys_solver.solve(rhs)

        # lift down vector
        self.uk_sys = temp[0:self.nv]

        self.k += 1
        self.t += self.dt

    def save(self):
        if (self.k-1) % const.LINEARIZED_SIM_SAVE_FREQ == 0:

            # compute norm of u_delta
            print self.t, np.linalg.norm(self.uk_sys)


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
            self.save()
            self.uk_next()

