from dolfin import *

class Newton():
    """Newton Method for solving stationary navier stokes"""

    def __init__(self, const, ref, RE, REinitial=None, shear=1):
        # parameters
        self.const = const
        self.ref = ref
        self.RE = RE
        self.shear = shear
        self.nu = const.STATIONARY_NU(RE)
        self.REinitial = REinitial

        # mesh and function spaces
        self.mesh = Mesh(self.const.MESH_XML(ref))
        self.V = VectorFunctionSpace(self.mesh, self.const.STATIONARY_V, self.const.STATIONARY_V_DIM)
        self.Q = FunctionSpace(self.mesh, self.const.STATIONARY_Q, self.const.STATIONARY_Q_DIM)
        self.W = MixedFunctionSpace([self.V, self.Q])

        # right hand side (external force)
        self.rhs = self.const.STATIONARY_RHS

    def solve(self):

        # define test functions
        (v, q) = TestFunctions(self.W)

        # define trial function
        if self.REinitial:
            w = Function(self.W, self.const.STATIONARY_W_XML(self.ref, self.REinitial))
        else:
            w = Function(self.W)

        (u, p) = (as_vector((w[0], w[1])), w[2])

        # get boundary conditions
        bc = self.const.STATIONARY_BOUNDARY_CONDITIONS(self.W)

        # build weak formulation
        a1 = inner(grad(u) * u, v) * dx
        a2 = self.nu * inner(grad(u), grad(v)) * dx
        a3 = -1 * p * div(v) * dx
        cond = -1 * div(u) * q * dx
        rhs = inner(self.rhs, v) * dx
        F = a1 + a2 + a3 + cond + rhs
        n = FacetNormal(self.mesh)


        # build derivative, problem and solver
        dw = TrialFunction(self.W)
        dF = derivative(F, w, dw)
        nsproblem = NonlinearVariationalProblem(F, w, bc, dF)



        # add optimality constraint
        if self.shear == 0:
            solver = NonlinearVariationalSolver(nsproblem)
            solver.parameters["newton_solver"]["maximum_iterations"] = self.const.STATIONARY_NEWTON_STEPS
            solver.parameters["newton_solver"]["absolute_tolerance"] = self.const.STATIONARY_NEWTON_ABS_TOL
            solver.parameters["newton_solver"]["relative_tolerance"] = self.const.STATIONARY_NEWTON_REL_TOL
            solver.solve()

        else:
            if self.shear == 1:
                shearfunction = MeshFunction("size_t", self.mesh, self.const.SHEAR_XML(self.ref, 1))
                ds = Measure("ds")[shearfunction]
                M = inner(dot(grad(u), Constant((0.0, -1.0))), dot(grad(u), Constant((0.0, -1.0))))*ds(self.const.GAMMASHEAR1_INDICES)
            elif self.shear == 2:
                shearfunction = MeshFunction("size_t", self.mesh, self.const.SHEAR_XML(self.ref, 2))
                ds = Measure("ds")[shearfunction]
                M = inner(dot(grad(u), Constant((0.0, -1.0))), dot(grad(u), Constant((0.0, -1.0))))*ds(self.const.GAMMASHEAR2_INDICES)
            elif self.shear == 3:
                shearfunction = MeshFunction("size_t", self.mesh, self.const.SHEAR_XML(self.ref, 3))
                ds = Measure("ds")[shearfunction]
                M = inner(dot(grad(u), Constant((0.0, -1.0))), dot(grad(u), Constant((0.0, -1.0))))*ds(self.const.GAMMASHEAR3_INDICES)


            solver = AdaptiveNonlinearVariationalSolver(nsproblem, M)
            solver.parameters["max_iterations"] = self.const.STATIONARY_ADAPTIVE_STEPS
            solver.parameters["nonlinear_variational_solver"]["newton_solver"]["maximum_iterations"] = self.const.STATIONARY_NEWTON_STEPS
            solver.parameters["nonlinear_variational_solver"]["newton_solver"]["absolute_tolerance"] = self.const.STATIONARY_NEWTON_ABS_TOL
            solver.parameters["nonlinear_variational_solver"]["newton_solver"]["relative_tolerance"] = self.const.STATIONARY_NEWTON_REL_TOL
            solver.solve(self.const.STATIONARY_ADAPTIVE_TOL)


        # split w
        if self.shear==0:
            (u, p) = w.split(deepcopy=True)
        else:
            (u, p) = w.leaf_node().split(deepcopy=True)
        self.u = interpolate(u, self.V)
        self.p = interpolate(p, self.Q)
        self.w = w

        # Compute divergence
        #self._check_divergence()
        #self._checku2()
        #self._checkdiff()


    def _check_divergence(self):
        """Check divergence of velocity."""

        # Compute L2 norm of divergence
        print "||div u||_L2 =", norm(self.u, "Hdiv0")

        # Compute projection of div u into Q_0
        pdivu = project(div(self.u), self.Q)
        zero = Constant(0.0)
        bc = DirichletBC(self.Q, zero, DomainBoundary())
        bc.apply(pdivu.vector())

        # Compute "weak" L2 norm of divergence
        print "||div u||_w  =", sqrt(abs(assemble(pdivu * div(self.u) * dx)))

    def _checku2(self):
        print "||u2||_L2=", norm(self.u2, "L2")
        print "||u2||_H1=", norm(self.u2, "H1")
        print "||u2||_H10=",  norm(self.u2, "H10")
        print "||u2||_Hdiv=", norm(self.u2, "Hdiv")
        print "||u2||_Hdiv0=", norm(self.u2, "Hdiv0")

    def _checkdiff(self):
        print "||u2||_L2=", errornorm(self.u, self.u2, "L2")
        print "||u2||_H1=", errornorm(self.u, self.u2, "H1")
        print "||u2||_H10=", errornorm(self.u, self.u2,  "H10")
        print "||u2||_Hdiv=", errornorm(self.u, self.u2, "Hdiv")
        print "||u2||_Hdiv0=", errornorm(self.u, self.u2, "Hdiv0")


    def save(self):
        File(self.const.STATIONARY_U_PVD(self.ref, self.RE)) << self.u
        File(self.const.STATIONARY_U_XML(self.ref, self.RE)) << self.u
        File(self.const.STATIONARY_P_PVD(self.ref, self.RE)) << self.p
        File(self.const.STATIONARY_P_XML(self.ref, self.RE)) << self.p
        File(self.const.STATIONARY_W_XML(self.ref, self.RE)) << self.w
