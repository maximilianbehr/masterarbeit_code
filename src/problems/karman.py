# -*- coding: utf-8 -*-
from problembase import ProblemBase
from problem_mesh.karman import GammaLower
from problem_mesh.karman import GammaLeft
from problem_mesh.karman import GammaUpper
from problem_mesh.karman import GammaBall
from problem_mesh.karman import GammaRight
from problem_mesh.karman import GammaBallCtrlLower
from problem_mesh.karman import GammaBallCtrlUpper
from dolfin import *


class Problem(ProblemBase):
    def __init__(self, options):
        ProblemBase.__init__(self, options)

        # load mesh
        self.mesh = Mesh(self.options["mesh"])

        # create right-hand side function
        self.f = Constant((0, 0))

        # set viscosity
        self.nu = 1.0 / options["RE"]

        # set end time
        self.T = options.get("T", None)


    def initial_conditions(V, Q):
        u0 = Constant((0, 0))
        p0 = Constant(0)
        return u0, p0

    def boundary_conditions(self, V, Q, t):
        # inflow condition
        self.u_in = Expression(("(1-x[1])*x[1]*(t/(1.0+t))", "0.0"), t=t)
        self.bcs_inflow = DirichletBC(V, self.u_in, GammaLeft())

        # no slip condition
        self.noslip = Constant((0.0, 0.0))
        self.bcs_upper = DirichletBC(V, self.noslip, GammaUpper())
        self.bcs_lower = DirichletBC(V, self.noslip, GammaLower())
        self.bcs_ball = DirichletBC(V, self.noslip, GammaBall())
        self.bcs_ctrllower = DirichletBC(V, self.noslip, GammaBallCtrlLower())
        self.bcs_ctrlupper = DirichletBC(V, self.noslip, GammaBallCtrlUpper())

        # collect boundary conditions
        bcu = [self.bcs_upper, self.bcs_lower, self.bcs_ball, self.bcs_inflow, self.bcs_ctrllower, self.bcs_ctrlupper]

        # collect boundary conditions for pressure
        bcp = []

        return bcu, bcp

    def stat_boundary_conditions(self, W):

        # inflow condition
        self.u_in = Expression(("(1-x[1])*x[1]", "0.0"))
        self.bcs_inflow = DirichletBC(W.sub(0), self.u_in, GammaLeft())

        # no slip condition
        self.noslip = Constant((0.0, 0.0))
        self.bcs_upper = DirichletBC(W.sub(0), self.noslip, GammaUpper())
        self.bcs_lower = DirichletBC(W.sub(0), self.noslip, GammaLower())
        self.bcs_ball = DirichletBC(W.sub(0), self.noslip, GammaBall())
        self.bcs_ctrllower = DirichletBC(W.sub(0), self.noslip, GammaBallCtrlLower())
        self.bcs_ctrlupper = DirichletBC(W.sub(0), self.noslip, GammaBallCtrlUpper())

        # collect boundary conditions
        bcu = [self.bcs_upper, self.bcs_lower, self.bcs_ball, self.bcs_inflow, self.bcs_ctrllower, self.bcs_ctrlupper]

        # collect boundary conditions for pressure
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
