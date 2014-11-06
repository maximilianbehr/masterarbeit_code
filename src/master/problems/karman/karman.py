from ..problembase import *
from numpy import array
from mesh.mesh   import *

# Problem definition
class Problem(ProblemBase):

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        # Load mesh
        refinement_level = options["refinement_level"]
        if refinement_level > 5:
            raise RuntimeError, "No mesh available for refinement level %d" % refinement_level

        self.mesh = Mesh(KARMAN_REFN_MESH_FILE(parameters["refinement_algorithm"], refinement_level))


        # Create right-hand side function
        self.f =  Constant((0, 0))

        # Set viscosity (Re = 1000)
        self.nu = 1.0 / 1000.0

        # Characteristic velocity in the domain (used to determinde timestep)
        self.U = 3.5

        # Set end time
        self.T  = 1.0

    def initial_conditions(self, V, Q):

        u0 = Constant((0, 0))
        p0 = Constant(0)

        return u0, p0

    def boundary_conditions(self, V, Q, t):

        # Create boundary condition
        self.u_in               = Expression(("t*(1-x[1])*x[1]*2", "0.0"),t=t)
        self.u_inflow           = DirichletBC(V, self.u_in, GammaLeft())

        self.u_noslip           = Constant((0,0))
        self.u_noslip_upper     = DirichletBC(V, self.u_noslip, GammaUpper())
        self.u_noslip_lower     = DirichletBC(V, self.u_noslip, GammaLower())
        self.u_noslip_ball      = DirichletBC(V, self.u_noslip, GammaBall())

        # Collect boundary conditions
        bcu = [self.u_noslip_upper, self.u_noslip_lower, self.u_noslip_ball, self.u_inflow]
        bcp = []

        return bcu, bcp

    def update(self, t, u, p):
        self.u_in.t = t

    def functional(self, t, u, p):
        return 0.0

    def reference(self, t):
        return 0.0

    def __str__(self):
        return "Karman"
