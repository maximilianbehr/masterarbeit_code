from src.solvers.solverbase import *


class Solver(SolverBase):
    """Newton Method for solving stationary navier stokes"""

    def __init__(self, options):
        SolverBase.__init__(self, options)

    def solve(self, problem):

        # Get problem parameters
        mesh = problem.mesh

        # Define function spaces (P2-P1)
        self.V = VectorFunctionSpace(mesh, "CG", 2)
        self.Q = FunctionSpace(mesh, "CG", 1)
        W = self.V * self.Q

        # define test functions
        (v, q) = TestFunctions(W)

        # define trial function
        w = Function(W)
        (u, p) = (as_vector((w[0], w[1])), w[2])

        # rhs
        f = problem.f

        # Functions
        nu = Constant(problem.nu)

        # get boundary conditions
        bc = problem.stat_boundary_conditions(W)

        # build weak formulation
        a1 = inner(grad(u) * u, v) * dx
        a2 = nu * inner(grad(u), grad(v)) * dx
        a3 = -1 * p * div(v) * dx
        cond = -1 * div(u) * q * dx
        rhs = inner(f, v) * dx
        F = a1 + a2 + a3 + cond + rhs

        # build derivative
        dw = TrialFunction(W)
        dF = derivative(F, w, dw)

        # solve the problem
        nsproblem = NonlinearVariationalProblem(F, w, bc, dF)
        solver = NonlinearVariationalSolver(nsproblem)
        solver.parameters["newton_solver"]["maximum_iterations"] = self.options["newton_solver_max_iterations"]
        solver.parameters["newton_solver"]["absolute_tolerance"] = self.options["newton_solver_absolute_tolerance"]
        solver.parameters["newton_solver"]["relative_tolerance"] = self.options["newton_solver_relative_tolerance"]
        solver.solve()

        # split w
        (u, p) = w.split(deepcopy=True)

        # save
        self.save(u, p)

        # Compute divergence
        if self.options["compute_divergence"]:
            check_divergence(u, p.function_space())

        return u, p


    def eval(self):
        return 0, 0


    def save(self, u, p):
        interu = interpolate(u, self.V)
        interp = interpolate(p, self.Q)

        if self.options["u_pvd"]:
            File(self.options["u_pvd"]) << interu

        if self.options["p_pvd"]:
            File(self.options["p_pvd"]) << interp

        if self.options["u_xml"]:
            File(self.options["u_xml"]) << interu

        if self.options["p_xml"]:
            File(self.options["p_xml"]) << interp


    def __str__(self):
        return "Newton Stationary"
