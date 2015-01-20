import src.karman_const as const
from src.mesh.karman import GammaBallCtrlUpper
from src.mesh.karman import GammaBallCtrlLower
from src.mesh.karman import GammaBall
from src.mesh.karman import GammaLeft
from src.mesh.karman import GammaRight
from src.mesh.karman import GammaLower
from src.mesh.karman import GammaUpper
import scipy.io as scio
import scipy.sparse as scsp
import scipy.sparse.linalg as scspli
import numpy as np

from dolfin import *


class Linearized():
    def __init__(self, ref, RE, pertubationeps, dt, T):

        # parameters
        self.ref = ref
        self.RE = RE
        self.dt = dt
        self.T = T
        self.pertubationeps = pertubationeps
        self.t = 0
        self.k = 1

        # mesh and function spaces
        self.mesh = Mesh(const.MESH_XML(ref))
        self.V = VectorFunctionSpace(self.mesh, const.LINEARIZED_V, const.LINEARIZED_V_DIM)
        self.Q = FunctionSpace(self.mesh, const.LINEARIZED_Q, const.LINEARIZED_Q_DIM)
        self.u_stat = Function(self.V, const.STATIONARY_U_XML(ref, RE))

        # define attributes v_inner_nodes and v_outer_nodes
        self._inner_outer_nodes()

        # define state
        self.u_k_sys = self.u_stat.vector().array().reshape(len(self.u_stat.vector().array()), )[self.v_inner_nodes]
        self.u_k_sys *= self.pertubationeps

        # read matrices
        M = scio.mmread(const.LQR_M_MTX(ref, RE)).tocsr()
        S = scio.mmread(const.LQR_S_MTX(ref, RE)).tocsr()
        R = scio.mmread(const.LQR_R_MTX(ref, RE)).tocsr()
        K = scio.mmread(const.LQR_K_MTX(ref, RE)).tocsr()
        G = scio.mmread(const.LQR_G_MTX(ref, RE)).tocsr()
        Mlower = scio.mmread(const.LQR_MLOWER_MTX(ref, RE)).tocsr()
        Mupper = scio.mmread(const.LQR_MUPPER_MTX(ref, RE)).tocsr()
        B = scio.mmread(const.LQR_B_MTX(ref, RE))

        # compress system
        rowcps = lambda mat: mat.tocsr()[self.v_inner_nodes, :]
        colcps = lambda mat: mat.tocsc()[:, self.v_inner_nodes]
        self.Mcps = colcps(rowcps(M)).tocsc()
        self.Scps = colcps(rowcps(S)).tocsc()
        self.Rcps = colcps(rowcps(R)).tocsc()
        self.Kcps = colcps(rowcps(K)).tocsc()
        self.Mlowercps = colcps(rowcps(Mlower)).tocsc()
        self.Muppercps = colcps(rowcps(Mupper)).tocsc()
        self.Gcps = rowcps(G).tocsc()
        self.Bcps = B[self.v_inner_nodes, :]

        # system sizes
        self.nv = self.Mcps.shape[0]
        self.np = self.Gcps.shape[1]

        # build system matrices
        self.Asys = -self.Scps-self.Rcps-self.Kcps
        #self.Asys -= self.Mlowercps-self.Muppercps
        self.Bsys = np.vstack((self.Bcps, np.zeros((self.np, 2))))
        u = scsp.hstack([self.Mcps - self.dt*self.Asys, self.dt*(-self.Gcps)])
        l = scsp.hstack([self.dt * (-self.Gcps.T), scsp.csr_matrix((self.np, self.np))])
        self.Msys_ode = scsp.vstack([u, l]).tocsc()
        self.Msys_lift = scsp.vstack([self.Mcps, scsp.csr_matrix((self.np, self.nv))])

        # build solver
        self.Mcompress_solver = scspli.spilu(self.Msys_ode)

        # visualisation
        self.u_dolfin = Function(self.V)
        self.ufile = File(const.LINEARIZED_U_PVD(ref, RE))

    def _inner_outer_nodes(self):
        """compute inner and outer nodes, add the attribute v_inner_nodes and v_outer_nodes"""

        # build and collect boundary conditions
        noslip = Constant((0.0, 0.0))
        bcs_left = DirichletBC(self.V, noslip, GammaLeft())
        bcs_upper = DirichletBC(self.V, noslip, GammaUpper())
        bcs_lower = DirichletBC(self.V, noslip, GammaLower())
        bcs_ball = DirichletBC(self.V, noslip, GammaBall())
        bcs_ctrllower = DirichletBC(self.V, noslip, GammaBallCtrlLower())
        bcs_ctrlupper = DirichletBC(self.V, noslip, GammaBallCtrlUpper())
        bcu = [bcs_left, bcs_upper, bcs_lower, bcs_ball, bcs_ctrllower, bcs_ctrlupper]

        # collect indices of outer nodes
        v_outer_nodes = []
        for bc in bcu:
            v_outer_nodes += bc.get_boundary_values().keys()
        v_outer_nodes = list(set(v_outer_nodes))

        # turn control on set node on the control as an inner node
        #self.v_outer_nodes = list(set(self.v_outer_nodes)-set(bcs_ctrllower.get_boundary_values().keys()))
        #self.v_outer_nodes = list(set(self.v_outer_nodes)-set(bcs_ctrlupper.get_boundary_values().keys()))

        # collect indices of inner nodes
        v_inner_nodes = list(set(range(0, self.V.dim()))-set(v_outer_nodes))

        self.v_inner_nodes = v_inner_nodes
        self.v_outer_nodes = v_outer_nodes

    def u_k_uncompress(self):
        """uncompress a state to all nodes, add zeros at the outer nodes"""
        u_k = np.zeros((self.V.dim(),))
        u_k[self.v_inner_nodes] = self.u_k_sys
        return u_k

    def assembleN(self):

        # Function and TestFunction
        w_test = TestFunction(self.V)

        # insert values at inner nodes
        u_k = self.u_k_uncompress()

        # assemble vector
        self.u_dolfin.vector().set_local(u_k)
        N = inner(grad(self.u_dolfin) * self.u_dolfin, w_test)*dx
        Nassemble = assemble(N)
        self.Nukuk = Nassemble.array().reshape((self.V.dim(),))

        # compress
        self.Nukuk = self.Nukuk[self.v_inner_nodes]

        # lift
        self.Nukuk = np.append(self.Nukuk, np.zeros((self.np,)))


    def u_k_next(self):

        # compute new rhs
        self.assembleN()
        rhs = self.Msys_lift * self.u_k_sys - self.dt*self.Nukuk

        #turn control on
        #Rhs -= (self.Bcompress[:,0]*0.01*sin(self.t)+ self.Bcompress[:,1]*0.02*cos(self.t))

        # solve for the next step
        u_k_next = self.Mcompress_solver.solve(rhs)

        # lift down vector
        self.u_k_sys = u_k_next[0:self.nv]

        self.k += 1
        self.t += self.dt

    def save(self):
        if (self.k-1) % const.LINEARIZED_SAVE_FREQ == 0:

            # compute norm of u_delta
            print self.t, np.linalg.norm(self.u_k_sys)

            # uncompress and add stationary solution
            u_k_uncps = self.u_k_uncompress()
            v = u_k_uncps + self.u_stat.vector().array()

            # set it to a fenics function
            self.u_dolfin.vector().set_local(v)

            # save visualization
            self.ufile << self.u_dolfin


    def solve_ode(self):
        while self.t < (self.T + DOLFIN_EPS):
            self.save()
            self.u_k_next()
