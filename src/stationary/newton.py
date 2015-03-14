from dolfin import *



class Newton():
    """Newton Method for solving stationary navier stokes"""

    def __init__(self, const, ref, RE, initialw=None):
        # parameters
        self.const = const
        self.ref = ref
        self.RE = RE
        self.nu = const.GET_NU(RE)
        self.initialw = initialw

        # mesh and function spaces
        self.mesh = Mesh(self.const.MESH_XML(ref))
        self.V = VectorFunctionSpace(self.mesh, self.const.V, self.const.V_DIM)
        self.Q = FunctionSpace(self.mesh, self.const.Q, self.const.Q_DIM)
        self.W = self.V*self.Q
        self.boundaryfunction = MeshFunction("size_t", self.mesh, const.BOUNDARY_XML(ref))

        # right hand side (external force)
        self.rhs = self.const.STATIONARY_RHS

    def solve(self):

        # define test functions
        (v, q) = TestFunctions(self.W)

        # define trial function
        if self.initialw:
            # w = Function(self.W, self.const.STATIONARY_W_XML(self.ref, self.REinitial))
            w = Function(self.W, self.initialw)
        else:
            w = Function(self.W)

        (u, p) = (as_vector((w[0], w[1])), w[2])

        # build weak formulation
        a1 = inner(grad(u) * u, v) * dx
        a2 = self.nu * inner(grad(u), grad(v)) * dx
        a3 = -1 * p * div(v) * dx
        cond = -1 * div(u) * q * dx
        rhs = inner(self.rhs, v) * dx
        F = a1 + a2 + a3 + cond + rhs

        # get boundary conditions
        bc = self.const.STATIONARY_BOUNDARY_CONDITIONS(self.W, self.boundaryfunction, self.const)

        # build derivative
        dw = TrialFunction(self.W)
        dF = derivative(F, w, dw)

        # solve the problem
        nsproblem = NonlinearVariationalProblem(F, w, bc, dF)
        solver = NonlinearVariationalSolver(nsproblem)

        solver.parameters["newton_solver"]["maximum_iterations"] = self.const.STATIONARY_NEWTON_STEPS
        solver.parameters["newton_solver"]["absolute_tolerance"] = self.const.STATIONARY_NEWTON_ABS_TOL
        solver.parameters["newton_solver"]["relative_tolerance"] = self.const.STATIONARY_NEWTON_REL_TOL
        solver.parameters["newton_solver"]["report"] = True
        solver.solve()

        # split w
        (u, p) = w.split(deepcopy=True)

        self.u = interpolate(u, self.V)
        self.p = interpolate(p, self.Q)
        self.w = w

        # Compute divergence
        self._check_divergence()

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


    def save(self):
        # set names and labels
        self.u.rename(*self.const.PVD_U_LABEL_NAME)
        self.p.rename(*self.const.PVD_P_LABEL_NAME)
        self.w.rename(*self.const.PVD_W_LABEL_NAME)

        # save pvd and xml files
        File(self.const.STATIONARY_U_PVD(self.ref, self.RE)) << self.u
        File(self.const.STATIONARY_U_XML(self.ref, self.RE)) << self.u
        File(self.const.STATIONARY_P_PVD(self.ref, self.RE)) << self.p
        File(self.const.STATIONARY_P_XML(self.ref, self.RE)) << self.p
        File(self.const.STATIONARY_W_XML(self.ref, self.RE)) << self.w
