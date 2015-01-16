from scipy.io import mmread
from scipy.sparse import hstack
from scipy.sparse import vstack
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve, splu, spilu
from dolfin import *
import numpy

from src.problems.problem_mesh.karman import GammaBall
from src.problems.problem_mesh.karman import GammaBallCtrlLower
from src.problems.problem_mesh.karman import GammaBallCtrlUpper
from src.problems.problem_mesh.karman import GammaLeft
from src.problems.problem_mesh.karman import GammaLower
from src.problems.problem_mesh.karman import GammaUpper
from src.problems.problem_mesh.karman import GammaRight

from dolfin import parameters
parameters["reorder_dofs_serial"] = False


class Linearized_NSE_SIM():
    def __init__(self, inoptions):

        # copy options
        self.options = inoptions.copy()

        # read mesh and build Function Spaces
        self.mesh = Mesh(self.options["mesh"])
        self.V = VectorFunctionSpace(self.mesh, "CG", 2)
        self.Q = FunctionSpace(self.mesh, "CG", 1)

        # build and collect boundary conditions
        noslip = Constant((0.0, 0.0))
        bcs_left = DirichletBC(self.V, noslip, GammaLeft())
        bcs_upper = DirichletBC(self.V, noslip, GammaUpper())
        bcs_lower = DirichletBC(self.V, noslip, GammaLower())
        bcs_ball = DirichletBC(self.V, noslip, GammaBall())
        bcs_ctrllower = DirichletBC(self.V, noslip, GammaBallCtrlLower())
        bcs_ctrlupper = DirichletBC(self.V, noslip, GammaBallCtrlUpper())
        self.bcu = [bcs_left, bcs_upper, bcs_lower, bcs_ball, bcs_ctrllower, bcs_ctrlupper]


        # collect indices of outer nodes
        self.v_outer_nodes = []
        for bc in self.bcu:
            self.v_outer_nodes += bc.get_boundary_values().keys()
        self.v_outer_nodes = list(set(self.v_outer_nodes))

        # turn control on set node on the control as an inner node
        self.v_outer_nodes = list(set(self.v_outer_nodes)-set(bcs_ctrllower.get_boundary_values().keys()))
        self.v_outer_nodes = list(set(self.v_outer_nodes)-set(bcs_ctrlupper.get_boundary_values().keys()))

        #import ipdb
        #ipdb.set_trace()


        # collect indices of inner nodes
        self.v_inner_nodes = list(set(range(0, self.V.dim()))-set(self.v_outer_nodes))

        # read matrices
        M = mmread(self.options["M_mtx"]).tocoo()
        S = mmread(self.options["S_mtx"]).tocoo()
        K = mmread(self.options["K_mtx"]).tocoo()
        R = mmread(self.options["R_mtx"]).tocoo()
        G = mmread(self.options["G_mtx"]).tocoo()
        Gt = mmread(self.options["Gt_mtx"]).tocoo()

        # turn control on
        B = mmread(self.options["B_mtx"])
        Mlower = mmread(self.options["Mlower_mtx"]).tocoo()
        Mupper = mmread(self.options["Mupper_mtx"]).tocoo()


        # set parameters
        self.dt = self.options["dt"]
        self.T = self.options["T"]
        self.t = 0
        self.k = 1

        # compress system
        self.Mcompress = M.tocsr()[self.v_inner_nodes,:].tocsc()[:,self.v_inner_nodes].tocsc()

        #self.Acompress = (-S - R - K).tocsr()[self.v_inner_nodes,:].tocsc()[:,self.v_inner_nodes].tocsc()
        # turn control on
        self.Acompress = (-S-R-K-Mlower-Mupper).tocsr()[self.v_inner_nodes,:].tocsc()[:,self.v_inner_nodes].tocsc()
        self.Gcompress = G.tocsr()[self.v_inner_nodes,:].tocsc()

        self.nv = self.Mcompress.shape[0]
        self.np = self.Gcompress.shape[1]

        # turon control on
        Bcompress = B[self.v_inner_nodes,:]
        self.Bcompress = numpy.vstack((Bcompress, numpy.zeros((self.np,2))))

        #import ipdb
        #ipdb.set_trace()

        # build block matrices
        upperblock = hstack([self.Mcompress - self.dt*self.Acompress, self.dt*(-self.Gcompress)])
        lowerblock = hstack([self.dt * (-self.Gcompress.T), csr_matrix((self.np, self.np))])
        self.Mcompress_ode = vstack([upperblock, lowerblock])
        self.Mcompress_lift = vstack([self.Mcompress, csr_matrix((self.np, self.nv))])

        # build solver
        self.Mcompress_solver = spilu(self.Mcompress_ode.tocsc())

        # read stationary solution
        self.u_stat = Function(self.V, self.options["u_stat"])

        # conversions
        self.u_k_compress = self.u_stat.vector().array().reshape(len(self.u_stat.vector().array()), )[self.v_inner_nodes]
        self.u_k_compress *= self.options["pertubation_eps"]

        # visualisation
        self.u_dolfin = Function(self.V)
        self._ufile = None

    def assemble_N(self):

        # Function and TestFunction
        w_test = TestFunction(self.V)

        # inser values at inner nodes
        u_k = self.u_k_uncompress()

        # assemble vector
        self.u_dolfin.vector().set_local(u_k)
        N = inner(grad(self.u_dolfin) * self.u_dolfin, w_test)*dx
        Nassemble = assemble(N)
        self.Nukuk = Nassemble.array().reshape((self.V.dim(),))

        # compress
        self.Nukuk = self.Nukuk[self.v_inner_nodes]

        # lift
        self.Nukuk = numpy.append(self.Nukuk, numpy.zeros((self.np,)))


    def u_k_next(self):

        # compute new rhs
        self.assemble_N()
        Rhs = self.Mcompress_lift * self.u_k_compress - self.dt*self.Nukuk

        #turn control on
        Rhs -= (self.Bcompress[:,0]*0.01*sin(self.t)+ self.Bcompress[:,1]*0.02*cos(self.t))

        # solve for the next step
        u_k_next = self.Mcompress_solver.solve(Rhs)

        # lift down vector
        self.u_k_compress = u_k_next[0:self.nv]

        self.k += 1
        self.t+=self.dt

    def u_k_uncompress(self):

        # add boundary conditions all zero
        u_k = numpy.zeros((self.V.dim(),))
        u_k [self.v_outer_nodes] = 0.0
        #for bc in self.bcu:
        #    bcd = bc.get_boundary_values()
        #    u_k[bcd.keys()]=bcd.values()

        u_k[self.v_inner_nodes] = self.u_k_compress

        return u_k

    def save(self):

        frequency = self.options["save_frequency"]
        if ((self.k - 1) % frequency == 0) or (self.options["save_solution_at_t=T"] and self.t >= self.T):
            print self.t, numpy.linalg.norm(self.u_k_compress)
            u_k_uncompress = self.u_k_uncompress()
            v = u_k_uncompress + self.u_stat.vector().array()

            self.u_dolfin.vector().set_local(v)


            if self.options["u_pvd"]:
                if self._ufile is None:
                    self._ufile = File(self.options["u_pvd"])
                self._ufile << self.u_dolfin

            if self.options["u_t_xml"]:
                file = File(self.options["u_t_xml"].format(self.t), "compressed")
                file << self.u_dolfin.vector()

    def solve_ode(self):
        while self.t < (self.T + DOLFIN_EPS):
            self.save()
            self.u_k_next()
