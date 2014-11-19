# -*- coding: utf-8 -*-

from src.solvers.solverbase import *


class Solver(SolverBase):
    "Consistent splitting scheme by Guermond and Shen."

    def __init__(self, options, order=2):
        SolverBase.__init__(self, options)
        self.order = order

    def solve(self, problem):

        # Get problem data
        mesh = problem.mesh
        dt, t, t_range = problem.timestep(problem)

        # Define function spaces
        V = VectorFunctionSpace(mesh, "CG", 2)
        Q = FunctionSpace(mesh, "CG", 1)

        # Get initial and boundary conditions
        u0, p0 = problem.initial_conditions(V, Q)
        bcu, bcp = problem.boundary_conditions(V, Q, t)
        bcpsi = homogenize(bcp)
        pbar = problem.pressure_bc(Q)

        # Remove boundary stress term is problem is periodic
        if is_periodic(bcp):
            beta = Constant(0)
        else:
            beta = Constant(1)

        # Test and trial functions
        v = TestFunction(V)
        q = TestFunction(Q)
        u = TrialFunction(V)
        p = TrialFunction(Q)

        # Functions
        u0 = interpolate(u0, V)
        u1 = interpolate(u0, V)
        p0 = interpolate(p0, Q)
        p1 = interpolate(p0, Q)
        p2 = interpolate(p0, Q)
        nu = Constant(problem.nu)
        k = Constant(dt)
        n = FacetNormal(mesh)
        f = problem.f
        psi = Function(Q)

        # Tentative pressure
        if self.order == 1:
            ps = p1
        else:
            ps = 2 * p1 - p0

        # Tentative velocity step
        F1 = (1 / k) * inner(v, u - u0) * dx + inner(v, grad(u0) * u0) * dx \
             + inner(epsilon(v), sigma(u, ps, nu)) * dx \
             - beta * nu * inner(v, grad(u).T * n) * ds + inner(v, pbar * n) * ds \
             - inner(v, f) * dx
        a1 = lhs(F1)
        L1 = rhs(F1)

        # Pressure correction
        a2 = inner(grad(q), grad(p)) * dx
        L2 = (1 / k) * inner(grad(q), u1 - u0) * dx - (1 / k) * inner(q * n, u1 - u0) * ds

        # Pressure update
        a3 = q * p * dx
        L3 = q * (ps + psi - nu * div(u1)) * dx

        # Assemble matrices
        A1 = assemble(a1)
        A2 = assemble(a2)
        A3 = assemble(a3)

        # Time loop
        prec = "amg" if has_krylov_solver_preconditioner("amg") else "default"

        self.start_timing()
        for t in t_range:

            # Get boundary conditions
            bcu, bcp = problem.boundary_conditions(V, Q, t)

            # Compute tentative velocity step
            b = assemble(L1)
            [bc.apply(A1, b) for bc in bcu]
            solve(A1, u1.vector(), b, "gmres", "ilu")

            # Compute pressure correction
            b = assemble(L2)
            if len(bcp) == 0 or is_periodic(bcp): normalize(b)
            [bc.apply(A2, b) for bc in bcpsi]
            if is_periodic(bcp):
                solve(A2, psi.vector(), b)
            else:
                solve(A2, psi.vector(), b, "gmres", prec)
            if len(bcp) == 0 or is_periodic(bcp): normalize(psi.vector())

            # Compute updated pressure
            b = assemble(L3)
            if len(bcp) == 0: normalize(b)
            [bc.apply(A3, b) for bc in bcp]
            solve(A3, p2.vector(), b, "gmres", "ilu")

            # Update
            self.update(problem, t, u1, p1)
            u0.assign(u1)
            p0.assign(p1)
            p1.assign(p2)

        return u1, p2

    def __str__(self):
        return "CSS"