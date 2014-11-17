from src.solvers.solverbase import *


class Solver(SolverBase):
    "Newton Method for solving stationary navier stokes"

    def __init__(self, options):
        SolverBase.__init__(self, options)


    def solve(self, problem):

        # Get problem parameters
        mesh = problem.mesh


        # Define function spaces (P2-P1)
        V = VectorFunctionSpace(mesh, "Lagrange", 2)
        Q = FunctionSpace(mesh, "Lagrange", 1)
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
        (u, p) = w.split()

        #save data
        if self.options["save_solution"]:
            # Save velocity and pressure
            solvername = self.__module__.split(".")[-1].lower()
            dir = problem.output_location(solvername)

            # Create files for saving
            print dir
            File(dir + "/velo.pvd") << u
            File(dir + "/pressure.pvd") << p

        # Save vectors in xml format
        if self.options["save_xml"]:
            file = File(dir + "/w_velo_pressure.xml", "compressed")
            file << w.vector()

        return Function(u), Function(p)


    def eval(self):
        return 0, 0

    def __str__(self):
        return "Newton Stationary"
