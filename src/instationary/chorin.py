# -*- coding: utf-8 -*-

from dolfin import *

import src.benchmarks.karman.karman_const as const
from aux import timestep
from aux import boundary_conditions
from aux import check_divergence


class Chorin():
    """Original pressure-correction scheme by Chorin and Temam."""

    def __init__(self, ref, RE, T, savefreq=None, dt=None):
        # parameters
        self.ref = ref
        self.RE = RE
        self.T = T
        self.dt = dt
        self.nu = const.INSTATIONARY_NU(RE)
        self.Umax = const.INSTATIONARY_UIN_MAX

        # mesh and function spaces
        self.mesh = Mesh(const.MESH_XML(ref))
        self.V = VectorFunctionSpace(self.mesh, const.INSTATIONARY_V, const.INSTATIONARY_V_DIM)
        self.Q = FunctionSpace(self.mesh, const.INSTATIONARY_Q, const.INSTATIONARY_Q_DIM)

        # right hand side (external forces)
        self.rhs = const.INSTATIONARY_RHS

        # files for saving
        self.savefreq = savefreq
        self.upvd = File(const.INSTATIONARY_U_PVD(ref, RE, "chorin"))
        self.ppvd = File(const.INSTATIONARY_P_PVD(ref, RE, "chorin"))

    def solve(self):

        #compute timesteps
        dt, t, t_range = timestep(self.T, self.Umax, self.nu, self.mesh, self.dt)

        # Get initial and boundary conditions
        u0 = const.INSTATIONARY_U0
        p0 = const.INSTATIONARY_P0

        bcu, bcp = boundary_conditions(self.V, self.Q, t)

        # Define test and trial functions
        v = TestFunction(self.V)
        q = TestFunction(self.Q)
        u = TrialFunction(self.V)
        p = TrialFunction(self.Q)

        # Define functions
        us = Function(self.V)
        u0 = interpolate(u0, self.V)
        u1 = interpolate(u0, self.V)
        p1 = interpolate(p0, self.Q)
        nu = Constant(self.nu)
        k = Constant(dt)
        f = self.rhs

        # Tentative velocity step
        F1 = (1.0 / k) * inner(v, u - u0) * dx + inner(v, grad(u0) * u0) * dx \
             + nu * inner(grad(v), grad(u)) * dx - inner(v, f) * dx
        a1 = lhs(F1)
        L1 = rhs(F1)

        # Poisson problem for the pressure
        a2 = inner(grad(q), grad(p)) * dx
        L2 = -(1 / k) * q * div(us) * dx

        # Velocity update
        a3 = inner(v, u) * dx
        L3 = inner(v, us) * dx - k * inner(v, grad(p1)) * dx

        # Assemble matrices
        A1 = assemble(a1)
        A2 = assemble(a2)
        A3 = assemble(a3)

        # Time loop
        prec = "amg" if has_krylov_solver_preconditioner("amg") else "default"
        step = 0
        for t in t_range:

            # Compute tentative velocity
            bcu, bcp = boundary_conditions(self.V, self.Q, t)
            b = assemble(L1)
            [bc.apply(A1, b) for bc in bcu]
            solve(A1, us.vector(), b, "gmres", "ilu")

            # Compute p1
            b = assemble(L2)
            if len(bcp) == 0:
                normalize(b)

            [bc.apply(A2, b) for bc in bcp]
            solve(A2, p1.vector(), b, "cg", prec)

            if len(bcp) == 0:
                normalize(p1.vector())

            # Compute u1
            b = assemble(L3)
            [bc.apply(A3, b) for bc in bcu]
            solve(A3, u1.vector(), b, "gmres", "ilu")

            # save
            if step % self.savefreq == 0:
                self.upvd << u1
                self.ppvd << p1

                print "{0:g}% done (t = {1:g}, T = {2:g})".format(100.0 * (t / self.T), t, self.T)

                check_divergence(u1, self.V, self.Q)


            u0.assign(u1)
            step+=1


    def __str__(self):
        return "Chorin"
