# -*- coding: utf-8 -*-


from numpy import linspace

from dolfin import *


class ProblemBase:
    "Base class for all problems."

    def __init__(self, options):

        # Store options
        self.options = options

        # Parameters must be defined by subclass
        self.mesh = None
        self.f = None
        self.bcu = []
        self.bcp = []
        self.nu = None
        self.t = 0
        self.T = None
        self.dt = None
        self.u0 = None
        self.p0 = None
        self.u = None
        self.p = None
        self.U = 1.0

        # FIXME: Temporary while testing
        self.bcpsi = []

    def output_location(self, solver):
        """output dir"""
        pass

    def update_problem(self, t, u, p):
        "Update problem at time t"

        # Update state
        self.t = t
        self.u = u
        self.p = p

        # Call problem-specific update
        self.update(t, u, p)

    def update(self, t, u, p):
        "Problem-specific update at time t"
        pass

    def functional(self, t, u, p):
        "Return value of functional of interest"
        return 0.0

    def reference(self, t):
        "Return reference value for functional"
        return None

    def tolerance(self, problem):
        "Return tolerance (used as local convergence criterion)."
        if str(problem) == 'Channel':
            return 1e-11
        elif str(problem) == 'Cylinder':
            return 1e-7
        else:
            return 1e-6

    # FIXME: Move to solverbase.py and rename select_timestep

    def timestep(self, problem):
        "Return time step and number of time steps for problem."

        # FIXME: This looks very complex, should be cleaned up

        T = self.options.get("T", problem.T)
        U = problem.U
        nu = problem.nu
        h = self.mesh.hmin()
        "Return time step and number of time steps for problem. Used for debugging / compilation only"
        if self.options["max_steps"] is not None:
            dt = 0.25 * h ** 2 / (U * (nu + h * U))
            n = self.options["max_steps"]
            T = n * dt
        else:

            if self.options["dt_division"] != 0 and self.options["dt"] is not None:
                dt = (self.options["dt"]) / int(sqrt(2) ** self.options["dt_division"])
                n = int(T / dt + 1.0)
                dt = T / n
                print 'Using options["dt"] and time step refinements'

            # Use time step specified in problem if available
            elif not self.options["dt"] is None:
                dt = self.options["dt"]
                print dt
                n = int(T / dt)
                print 'Using problem.dt'

            # Otherwise, base time step on problem_mesh size
            elif self.options["dt_division"] != 0:
                dt = 0.25 * h ** 2 / (U * (nu + h * U)) / int(sqrt(2) ** self.options["dt_division"])
                n = int(T / dt + 1.0)
                dt = T / n
                print 'Computing time step according to stability criteria and time step refinements'

            # Otherwise, base time step on problem_mesh size
            else:
                # dt =  0.25*h**2 / (U*(nu + h*U))
                dt = 0.2 * (h / U)
                n = int(T / dt + 1.0)
                dt = T / n
                print 'Computing time step according to stability criteria'

        # Compute range
        t_range = linspace(0, T, n + 1)[1:]  # FIXME: Comment out [1:] to run g2ref g2ref

        # Report time step
        print " "
        print 'Number of timesteps:', len(t_range)
        print 'Size of timestep:', dt
        print " "


        return dt, t_range[0], t_range

        # FIXME: Remove this and merge with boundary_condition()
        # def pressure_bc(self, Q):
        # return Constant(0)

