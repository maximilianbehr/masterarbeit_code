from src.solvers.solverbase import *


class Solver(SolverBase):
    "Newton Method for solving stationary navier stokes"

    def __init__(self, options):
        SolverBase.__init__(self, options)


    def solve(self, problem):

        # Get problem parameters
        mesh = problem.mesh


        # Define function spaces (P2-P1)
        V = VectorFunctionSpace(mesh, "CG", 2)
        Q = FunctionSpace(mesh, "CG", 1)
        W = V * Q

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
        solver.solve()

        # split w
        (u, p) = w.split(deepcopy=True)

        if self.options["sa"]:
            self.save(problem,u,p)


        return Function(u), Function(p)


    def eval(self):
        return 0, 0


           "RE": None,
           "u_pvd": None,
           "p_pvd": None,
           "u_xml": None,
           "p_xml": None,
           "debug": False,
           "krylov_solver_absolute_tolerance": 1e-25,
           "krylov_solver_relative_tolerance": 1e-12,
           "krylov_solver_monitor_convergence": False,
           }

    def save(self, problem, u, p):
        if self.options["u_pvd"]:
            File(self.options["u_pvd"])<<u

        if self.options["p_pvd"]:
            File(self.options["p_pvd"])<<p

        if self.options["u_xml"]:
            File(self.options["u_xml"])<<u

        if self.options["p_xml"]:
            File(self.options["p_xml"])<<p






    def __str__(self):
        return "Newton Stationary"
