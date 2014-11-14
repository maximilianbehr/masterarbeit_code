# -*- coding: utf-8 -*-

from problembase import *
import os

# Boundary value
class BoundaryValue(Expression):
    def value_shape(self):
        return (2,)
    def eval(self, values, x):
        if x[0] > DOLFIN_EPS and x[0] < 1.0 - DOLFIN_EPS and x[1] > 1.0 - DOLFIN_EPS:
            values[0] = 1.0
            values[1] = 0.0
        else:
            values[0] = 0.0
            values[1] = 0.0

# Problem definition
class Problem(ProblemBase):
    "2D lid-driven cavity test problem with known reference value."

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        # Create mesh
        N = options["N"]
        self.mesh = UnitSquare(N, N)

        # Create right-hand side function
        self.f = Constant((0, 0))

        # Set viscosity (Re = 1000)
        self.nu = 1.0 / 1000.0
        self.U = 1.0

        # Set end-time
        self.T = 2.5

    def output_location(self, solver):
        ref = options["refinement_level"]

        #build prefix of outputlocation
        EXTHARDDRIVE            = "/media/UNTITLED/"
        if socket.gethostname() == "pc747":
            raise NotImplementedError()
        elif socket.gethostname() == "pc800":
            if os.path.isdir(EXTHARDDRIVE):
                prefix = os.path.join(EXTHARDDRIVE,"results/karman")
            else:
                prefix = os.path.abspath("/scratch/behr/masters/src/master/results/karman/")
        elif socket.gethostname() == "pc785":
            raise NotImplementedError()
        elif socket.gethostname() == "pc633":
            raise NotImplementedError()
        elif socket.gethostname() == "jack":
            EXTHARDDRIVEMAC         = "/Volumes/UNTITLED/"
            if os.path.isdir(EXTHARDDRIVEMAC):
                prefix = os.path.join(EXTHARDDRIVEMAC,"data/karman")
            else:
                prefix = os.path.abspath("/Users/daniels/Documents/LiClipseWorkspace/master/src/master/data/karman")

        #build outputlocation
        return os.path.join(prefix,"%s/%s/RE_%.2e/%s/ref_%d"%(__version__,solver,self.nu,parameters["refinment_algorithm"],ref))



    def initial_conditions(self, V, Q):

        u0 = Constant((0, 0))
        p0 = Constant(0)

        return u0, p0

    def boundary_conditions(self, V, Q, t):

        element = VectorElement("CG", triangle, 1)
        self.g = BoundaryValue(element=element)
        bc = DirichletBC(V, self.g, DomainBoundary())

        return [bc], []

    def functional(self, t, u, p):

        # Only check final time
        if t < self.T:
            return 0
        else:
        # Compute stream function and report minimum
            psi = StreamFunction(u)
            vals  = psi.vector().array()
            vmin = vals.min()

            print "Stream function has minimal value" , vmin

            return vmin

    def reference(self, t):

        # Only check final time
        if t < self.T:
            return 0.0

        return -0.061076605

    def __str__(self):
        return "Driven cavity"

def StreamFunction(u):
    "Stream function for a given 2D velocity field."

    # Check dimension
    mesh = u.function_space().mesh()
    if not mesh.topology().dim() == 2:
        error("Stream-function can only be computed in 2D.")

    # Define variational problem
    V   = u.function_space().sub(0).collapse()
    q   = TestFunction(V)
    psi = TrialFunction(V)
    a   = dot(grad(q), grad(psi))*dx
    L   = dot(q, (u[1].dx(0) - u[0].dx(1)))*dx

    # Define boundary condition
    g  = Constant(0)
    bc = DirichletBC(V, g, DomainBoundary())

    # Compute solution
    psi = Function(V)
    solve(a == L, psi, bc)

    return psi
