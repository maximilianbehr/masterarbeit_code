# -*- coding: utf-8 -*-


from src.outputhandler import KarmanOutputHandler

from problembase import ProblemBase
from problem_mesh.karman import circle
from problem_mesh.karman import GammaLower
from problem_mesh.karman import GammaLeft
from problem_mesh.karman import GammaUpper
from problem_mesh.karman import GammaBall

from dolfin import *




# Problem definition
class Problem(ProblemBase):
    def __init__(self, options):
        ProblemBase.__init__(self, options)

        # Load problem_mesh
        self.mesh = Mesh(self.options["mesh"])

        # Create right-hand side function
        self.f = Constant((0, 0))

        # choose U such that U*2*r=1 and then RE=1/nu
        self.Umax = 1.0 / (2.0 * circle["r"])

        # Set viscosity (Re = 100)
        self.nu = 1.0/options["RE"]


    def RE(self):
        return (self.U * 2.0 *circle["r"]) / self.nu


    def initial_conditions(self, V, Q):

        u0 = Constant((0, 0))
        p0 = Constant(0)

        return u0, p0

    def boundary_conditions(self, V, Q, t):

        # Create boundary condition
        self.u_in = Expression(("U*(1-x[1])*x[1]*2*t", "0.0"),t=t , U=self.Umax)
        #self.u_in = Expression(('4*(x[1]*(1-x[1]))*sin(t*pi/8.0)', '0.0'),t=t)

        self.u_inflow = DirichletBC(V, self.u_in, GammaLeft())
        self.u_noslip = Constant((0, 0))
        self.u_noslip_upper = DirichletBC(V, self.u_noslip, GammaUpper())
        self.u_noslip_lower = DirichletBC(V, self.u_noslip, GammaLower())
        self.u_noslip_ball = DirichletBC(V, self.u_noslip,GammaBall())

        # Collect boundary conditions
        bcu = [self.u_noslip_upper, self.u_noslip_lower, self.u_noslip_ball, self.u_inflow]


        # boundary conditions for pressure
        # self.p_out              = Constant(0)
        # self.p_right            = DirichletBC(Q,self.p_out,GammaRight())
        # bcp = [self.p_right]
        bcp = []

        return bcu, bcp

    def stat_boundary_conditions(self, W):

        # Create boundary condition
        u_in = Expression(("U*(1-x[1])*x[1]*2", "0"), U=self.Umax)
        noslip_upper = DirichletBC(W.sub(0), (0, 0), GammaUpper())
        noslip_lower = DirichletBC(W.sub(0), (0, 0), GammaLower())
        noslip_ball = DirichletBC(W.sub(0), (0, 0), GammaBall())
        inflow = DirichletBC(W.sub(0), u_in, GammaLeft())
        bcu = [noslip_upper, noslip_lower, noslip_ball, inflow]


        # boundary conditions for pressure
        # self.p_out              = Constant(0)
        # self.p_right            = DirichletBC(Q,self.p_out,GammaRight())
        # bcp = [self.p_right]
        bcp = []

        return bcu + bcp


    def update(self, t, u, p):
        self.u_in.t = t

    def functional(self, t, u, p):
        return 0.0

    def reference(self, t):
        return 0.0

    def __str__(self):
        return "Karman"
