# -*- coding: utf-8 -*-


from problembase import *
from numpy import array

# Constants related to the geometry
bmarg = 1.e-3 + DOLFIN_EPS
xmin = 0.0
xmax = 2.2
ymin = 0.0
ymax = 0.41
xcenter = 0.2
ycenter = 0.2
radius = 0.05

# Inflow boundary
class InflowBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[0] < xmin + bmarg

# No-slip boundary
class NoslipBoundary(SubDomain):
    def inside(self, x, on_boundary):
        dx = x[0] - xcenter
        dy = x[1] - ycenter
        r = sqrt(dx*dx + dy*dy)
        return on_boundary and \
               (x[1] < ymin + bmarg or x[1] > ymax - bmarg or \
                r < radius + bmarg)

# Outflow boundary
class OutflowBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[0] > xmax - bmarg

# Problem definition
import os

class Problem(ProblemBase):

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        # Load mesh
        refinement_level = options["refinement_level"]
        if refinement_level > 5:
            raise RuntimeError, "No mesh available for refinement level %d" % refinement_level

        self.mesh = Mesh("data/cylinder_%d.xml.gz" % refinement_level)


        # Create right-hand side function
        self.f =  Constant((0, 0))

        # Set viscosity (Re = 1000)
        self.nu = 1.0 / 1000.0

        # Characteristic velocity in the domain (used to determinde timestep)
        self.U = 3.5

        # Set end time
        self.T  = 8.0

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

        # Create inflow boundary condition
        self.g0 = Expression(('4*Um*(x[1]*(ymax-x[1]))*sin(pi*t/8.0)/(ymax*ymax)', '0.0'),
                             Um=1.5, ymax=ymax, t=t)
        self.b0 = InflowBoundary()
        bc0 = DirichletBC(V, self.g0, self.b0)

        # Create no-slip boundary condition
        self.b1 = NoslipBoundary()
        self.g1 = Constant((0, 0))
        bc1     = DirichletBC(V, self.g1, self.b1)

        # Create outflow boundary condition for pressure
        self.b2 = OutflowBoundary()
        self.g2 = Constant(0)
        bc2     = DirichletBC(Q, self.g2, self.b2)

        # Collect boundary conditions
        bcu = [bc0, bc1]
        bcp = [bc2]

        return bcu, bcp

    def update(self, t, u, p):
        self.g0.t = t

    def functional(self, t, u, p):

        if t < self.T:
            return 0.0

        x1 = array((xcenter - radius - DOLFIN_EPS, ycenter))
        x2 = array((xcenter + radius + DOLFIN_EPS, ycenter))

        return p(x1) - p(x2)

    def reference(self, t):

        if t < self.T:
            return 0.0

        return -0.111444953719

    def __str__(self):
        return "Cylinder"
