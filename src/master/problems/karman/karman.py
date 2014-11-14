# -*- coding: utf-8 -*-

from ..problembase import *
from mesh  import *
from dolfin import __version__

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
        self.nu = 1.0 / 100.0

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

        # Create boundary condition
        self.u_in               = Expression(("(1-x[1])*x[1]*2", "0.0"),t=t)
        #self.u_in = Expression(('4*(x[1]*(1-x[1]))*sin(pi*t/8.0)', '0.0'), t=t)


        self.u_inflow           = DirichletBC(V, self.u_in, GammaLeft())

        self.u_noslip           = Constant((0,0))
        self.u_noslip_upper     = DirichletBC(V, self.u_noslip, GammaUpper())
        self.u_noslip_lower     = DirichletBC(V, self.u_noslip, GammaLower())
        self.u_noslip_ball      = DirichletBC(V, self.u_noslip, GammaBall())

        # Collect boundary conditions
        bcu = [self.u_noslip_upper, self.u_noslip_lower, self.u_noslip_ball, self.u_inflow]


        # boundary conditions for pressure
        #self.p_out              = Constant(0)
        #self.p_right            = DirichletBC(Q,self.p_out,GammaRight())
        #bcp = [self.p_right]
        bcp = []


        return bcu, bcp

    def stat_boundary_conditions(self, W):

        # Create boundary condition
        u_in                = Expression(("(1-x[1])*x[1]*2","0"))
        noslip_upper        = DirichletBC(W.sub(0), (0, 0), GammaUpper())
        noslip_lower        = DirichletBC(W.sub(0), (0, 0), GammaLower())
        noslip_ball         = DirichletBC(W.sub(0), (0, 0), GammaBall())
        inflow              = DirichletBC(W.sub(0),  u_in , GammaLeft())
        bcu                 = [noslip_upper,noslip_lower,noslip_ball,inflow]


        # boundary conditions for pressure
        #self.p_out              = Constant(0)
        #self.p_right            = DirichletBC(Q,self.p_out,GammaRight())
        #bcp = [self.p_right]
        bcp = []


        return bcu+bcp




    def update(self, t, u, p):
        self.u_in.t = t

    def functional(self, t, u, p):
        return 0.0

    def reference(self, t):
        return 0.0

    def __str__(self):
        return "Karman"
