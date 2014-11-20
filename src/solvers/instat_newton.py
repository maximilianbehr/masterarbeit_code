# -*- coding: utf-8 -*-

from src.solvers.solverbase import *
import warning
from dolfin.cpp.common import info

raise NotImplementedError("Code is not tested. There are faster solvers for solving instationary navier stokes. Dont waste my time for that.")


class Solver(SolverBase):
    """Newton Method plus theta time integration scheme for solving instationary navier stokes"""

    def __init__(self, options):
        SolverBase.__init__(self, options)

    def solve(self, problem):
        info("Newton Instationary assumes external forces f=0.")

        theta = 0.5 #crank nicelson

        # Get problem parameters
        mesh = problem.mesh

        if problem.dt is None:
            warning.warn("small steps sizes in time, solver needs long time,specifiy dt in your problem")

        warning.warn("solver assumes equidistant time discretization")
        dt, t, t_range = problem.timestep(problem)
        idt = 1.0/dt


        # Define function spaces (P2-P1)
        V = VectorFunctionSpace(mesh, "CG", 2)
        Q = FunctionSpace(mesh, "CG", 1)
        W = V * Q

         # define test functions
        (v, q) = TestFunctions(W)

        # define trial function
        w = Function(W)
        (u, p) = (as_vector((w[0], w[1])), w[2])

        # Get initial and boundary conditions
        bcu, bcp = problem.boundary_conditions(V, Q, t)
        bc = bcu + bcp


        # current master step
        T = sigma(u,p,problem.nu)


        # previous master step
        w0 = Function(W)
        (u0, p0) = (as_vector((w0[0], w0[1])), w0[2])
        T0 = sigma(u0,p,problem.nu)


        # Define variational forms without master derivative in previous master
        F0_eq1 = (inner(T0, grad(v))) * dx + inner(grad(u0) * u0, v) * dx
        F0_eq2 = 0 * q * dx
        F0 = F0_eq1 + F0_eq2

        # variational form without master derivative in current master
        F1_eq1 = (inner(T, grad(v)) + inner(grad(u) * u, v)) * dx
        F1_eq2 = q * div(u) * dx
        F1 = F1_eq1 + F1_eq2

        # combine variational forms with master derivative
        #
        # dw/dt + F(t) = 0 is approximated as
        # (w-w0)/dt + (1-theta)*F(t0) + theta*F(t) = 0
        #
        F = idt * inner((u - u0), v) * dx + (1.0 - theta) * F0 + theta * F1

        # residual of strong Navier-Stokes
        r = idt * (u - u0) + theta * grad(u) * u + (1.0 - theta) * grad(u0) * u0 \
            - theta * div(T) - (1.0 - theta) * div(T0)

        # define Jacobian
        J = derivative(F, w)

        # Time-stepping
        for t in t_range:

            # create variational problem and solver
            bcu, bcp = problem.boundary_conditions(V, Q, t)
            bc = bcu + bcp
            problem = NonlinearVariationalProblem(F, w, bc, J)
            solver = NonlinearVariationalSolver(problem)
            solver.parameters['newton_solver']['maximum_iterations'] = 20

            # Extract solutions:
            (u, p) = w.split()

            # update problem and save to file
            self.update(problem, t, u, p)

            # Move to next master step
            w0.assign(w)


        return u, p


    def __str__(self):
        return "Newton Instationary plus Theta Scheme"
