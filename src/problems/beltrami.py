# -*- coding: utf-8 -*-

from src.problems.problembase import *


# Problem definition
class Problem(ProblemBase):
    """3D test problem with known analytical solution."""

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        # We start with a UnitCube and modify it to get the problem_mesh we
        # want: (-1, 1) x (-1, 1) x (-1, 1)

        mesh_sizes = [5, 8, 11, 16, 23, 32]
        level = options["refinement_level"]
        N = int(mesh_sizes[level])

        self.mesh = UnitCubeMesh(N, N, N)
        self.scale = 2 * (self.mesh.coordinates() - 0.5)
        self.mesh.coordinates()[:, :] = self.scale

        # The body force term
        self.f = Constant((0, 0, 0))

        # Set the kinematic viscosity (nu = eta/rho)
        self.nu = 1.0

        # Set final time
        self.T = 0.5

        # FIXME: Write analytical solution in exact form as in paper

        # The analytical solution
        # Velocity
        self.analytical_u = \
            (
                '-((a*(pow(E,a*x[2])*cos(a*x[0] + d*x[1]) + pow(E,a*x[0])*sin(a*x[1] + d*x[2])))/pow(E,pow(d,2)*t*etabyrho))',
                '-((a*(pow(E,a*x[0])*cos(a*x[1] + d*x[2]) + pow(E,a*x[1])*sin(d*x[0] + a*x[2])))/pow(E,pow(d,2)*t*etabyrho))',
                '-((a*(pow(E,a*x[1])*cos(d*x[0] + a*x[2]) + pow(E,a*x[2])*sin(a*x[0] + d*x[1])))/pow(E,pow(d,2)*t*etabyrho))')
        # Pressure
        self.analytical_p = \
            (
                '-(rho/2.0)*(pow(a,2)*(pow(E,2*a*x[0]) + pow(E,2*a*x[1]) + pow(E,2*a*x[2]) + 2*pow(E,a*(x[1] + x[2]))*cos(d*x[0] + a*x[2])*sin(a*x[0] + d*x[1]) + 2*pow(E,a*(x[0] + x[1]))*cos(a*x[1] + d*x[2])*sin(d*x[0] + a*x[2]) + 2*pow(E,a*(x[0] + x[2]))*cos(a*x[0] + d*x[1])*sin(a*x[1] + d*x[2])))/(pow(E,pow(d,2)*t*etabyrho))')

        # Common parameters pertinent to the functional forms above
        self.u_params = {'a': pi / 4.0, 'd': pi / 2.0, 'E': e, 'etabyrho': 1.0, 't': 0.0}
        self.p_params = {'a': pi / 4.0, 'd': pi / 2.0, 'E': e, 'rho': 1.0, 'etabyrho': 1.0, 't': 0.0}


    def initial_conditions(self, V, Q):

        # Use analytical solutions at t = 0 as initial values
        self.exact_u = Expression(self.analytical_u, degree=3, **self.u_params)
        self.exact_p = Expression(self.analytical_p, degree=3, **self.p_params)

        return self.exact_u, self.exact_p

    def boundary_conditions(self, V, Q, t):
        self.exact_u.t = t
        self.exact_p.t = t

        bc0 = DirichletBC(V, self.exact_u, DomainBoundary())

        bcu = [bc0]
        bcp = []

        return bcu, bcp

    def update(self, t, u, p):
        print 'Time in update is:', t
        self.exact_u.t = t
        self.exact_p.t = t
        pass

    def functional(self, t, u, p):
        print 'Time in functional is:', t
        if t < self.T:
            return 0.0
        else:
            return errornorm(self.exact_u, u) / norm(self.exact_u, mesh=self.mesh)

    def reference(self, t):
        return 0.0

    def __str__(self):
        return "Beltrami"
