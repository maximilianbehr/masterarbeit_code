from scipy.io import mmread
from scipy.sparse import hstack
from scipy.sparse import vstack
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve
from dolfin import *
from src.aux import expand_vp_dolfunc


class Linearized_NSE_SIM():
    def __init__(self, inoptions):
        self.options = inoptions.copy()

        M = mmread(self.options["M_mtx"])
        S = mmread(self.options["S_mtx"])
        Mlower = mmread(self.options["Mlower_mtx"])
        Mupper = mmread(self.options["Mupper_mtx"])
        K = mmread(self.options["K_mtx"])
        R = mmread(self.options["R_mtx"])
        G = mmread(self.options["G_mtx"])
        Gt = mmread(self.options["Gt_mtx"])
        B = mmread(self.options["B_mtx"])

        self.np = G.shape[1]
        self.nv = M.shape[0]

        self.dt = self.options["dt"]
        self.T = self.options["T"]
        self.t = 0
        self.k = 1

        upperblock = hstack([M - self.dt * (-S - Mlower - Mupper - K - R), self.dt * (-G)])
        lowerblock = hstack([self.dt * (-Gt), csr_matrix((self.np, self.np))])
        self.M_ode = vstack([upperblock, lowerblock])
        self.M_lift = vstack([M, csr_matrix((self.np, self.nv))])

        self.mesh = Mesh(self.options["mesh"])
        self.V = VectorFunctionSpace(self.mesh, "CG", 2)
        self.Q = FunctionSpace(self.mesh, "CG", 1)

        self.u_stat = Function(self.V, self.options["u_stat"])

        import ipdb

        ipdb.set_trace()

        # conversions and scipy io functions
        self.u_k = self.u_stat.vector().array().reshape(len(self.u_stat.vector().vector()), 1)
        self.u_k *= self.options["pertubation_eps"]

    def u_k_next(self):
        Muk = self.M_lift * self.u_k
        u_k_next = spsolve(self.M_ode, Muk, permc_spec="COLAMD", use_umfpack="True")
        self.k += 1
        self.u_k = u_k_next[0:self.nv]

    def save(self):

        frequency = self.options["save_frequency"]
        if ((self.k - 1) % frequency == 0) or (self.options["save_solution_at_t=T"] and self.t >= self.T):
            # u to dolfin
            u_dolfin = expand_vp_dolfunc(V=self.V, Q=self.Q, invinds=None, diribcs=[], vp=None, vc=self.u_k, pc=None)

            if self.options["u_pvd"]:
                if self._ufile is None:
                    self._ufile = File(self.options["u_pvd"])
                self._ufile << (u_dolfin, self.t)

            if self.options["u_xml"]:
                file = File(self.options["u_t_xml"].format(self.t), "compressed")
                file << u_dolfin.vector()

    def solve_ode(self):
        while self.t < (self.T + DOLFIN_EPS):
            # self.save()
            self.u_k_next()
