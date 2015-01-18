import src.karman_const as const
from src.mesh.karman import GammaLeft
from src.mesh.karman import GammaLower
from src.mesh.karman import GammaUpper
from src.mesh.karman import GammaBall
from src.mesh.karman import GammaBallCtrlUpper
from src.mesh.karman import GammaBallCtrlLower

from dolfin import *



class Newton():
    """Newton Method for solving stationary navier stokes"""

    def __init__(self, ref, RE, REinitial=None):
        # parameters
        self.ref = ref
        self.RE = RE
        self.nu = const.STATIONARY_NU(RE)
        self.REinitial = REinitial

        # mesh and function spaces
        self.mesh = Mesh(const.MESH_XML(ref))
        self.V = VectorFunctionSpace(self.mesh, const.STATIONARY_V, const.STATIONARY_V_DIM)
        self.Q = FunctionSpace(self.mesh, const.STATIONARY_Q, const.STATIONARY_Q_DIM)
        self.W = self.V*self.Q


        # right hand side (external force)
        self.rhs = const.STATIONARY_RHS


    def _boundary_conditions(self):

        # inflow profile
        uin = const.STATIONARY_UIN

        # noslip at boundary parts
        noslip = Constant((0.0, 0.0))

        # define and collect boundary conditions
        bcu = [DirichletBC(self.W.sub(0), uin, GammaLeft()),
               DirichletBC(self.W.sub(0), noslip, GammaLower()),
               DirichletBC(self.W.sub(0), noslip, GammaUpper()),
               DirichletBC(self.W.sub(0), noslip, GammaBall()),
               DirichletBC(self.W.sub(0), noslip, GammaBallCtrlLower()),
               DirichletBC(self.W.sub(0), noslip, GammaBallCtrlUpper())]

        return bcu

    def solve(self):

        # define test functions
        (v, q) = TestFunctions(self.W)

        # define trial function
        if self.REinitial:
            w = Function(self.W, const.STATIONARY_W_XML(self.ref,self.REinitial))
        else:
            w = Function(self.W)

        (u, p) = (as_vector((w[0], w[1])), w[2])

        # get boundary conditions
        bc = self._boundary_conditions()

        # build weak formulation
        a1 = inner(grad(u) * u, v) * dx
        a2 = self.nu * inner(grad(u), grad(v)) * dx
        a3 = -1 * p * div(v) * dx
        cond = -1 * div(u) * q * dx
        rhs = inner(self.rhs, v) * dx
        F = a1 + a2 + a3 + cond + rhs

        # build derivative
        dw = TrialFunction(self.W)
        dF = derivative(F, w, dw)

        # solve the problem
        nsproblem = NonlinearVariationalProblem(F, w, bc, dF)
        solver = NonlinearVariationalSolver(nsproblem)

        solver.parameters["newton_solver"]["maximum_iterations"] = const.STATIONARY_NEWTON_STEPS
        solver.parameters["newton_solver"]["absolute_tolerance"] = const.STATIONARY_NEWTON_ABS_TOL
        solver.parameters["newton_solver"]["relative_tolerance"] = const.STATIONARY_NEWTON_REL_TOL
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
        File(const.STATIONARY_U_PVD(self.ref, self.RE)) << self.u
        File(const.STATIONARY_U_XML(self.ref, self.RE)) << self.u
        File(const.STATIONARY_P_PVD(self.ref, self.RE)) << self.p
        File(const.STATIONARY_P_XML(self.ref, self.RE)) << self.p
        File(const.STATIONARY_W_XML(self.ref, self.RE)) << self.w


    def __str__(self):
        return "Newton Stationary"
