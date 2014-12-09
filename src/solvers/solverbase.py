# -*- coding: utf-8 -*-
from os import getpid
from commands import getoutput

from dolfin import *


class SolverBase:
    "Base class for all solvers."

    def __init__(self, options):
        print "Dof reordering {0:s}".format(str(parameters["reorder_dofs_serial"]))

        # Store options
        self.options = options

        # Reset some solver variables
        self._time = None
        self._cputime = 0.0
        self._timestep = 0

        # Reset files for storing solution
        self._ufile = None
        self._pfile = None

        # Reset storage for functional values and errors
        self._t = []
        self._M = []
        self._m = []
        self._e = []

    def getMyMemoryUsage(self):
        mypid = getpid()
        mymemory = getoutput("ps -o rss %s" % mypid).split()[1]
        return mymemory

    def start_timing(self):
        """Start timing, will be paused automatically during update
        and stopped when the end-time is reached."""
        self._time = time()

    def solve(self, problem, dt, plot_solution=True):
        "Solve problem"
        raise NotImplementedError


    def save(self, problem, t, u, p):

        frequency = self.options["save_frequency"]
        if ((self._timestep - 1) % frequency == 0) or (self.options["save_solution_at_t=T"] and t >= problem.T):

            if self.options["u_pvd"]:
                if self._ufile is None:
                    self._ufile = File(self.options["u_pvd"])
                self._ufile << (u, t)

            if self.options["p_pvd"]:
                if self._pfile is None:
                    self._pfile = File(self.options["p_pvd"])
                self._pfile << (p, t)

            if self.options["u_xml"]:
                file = File(self.options["u_xml"].format(t), "compressed")
                file << u.vector()

            if self.options["p_xml"]:
                file = File(self.options["p_xml"].format(t), "compressed")
                file << p.vector()


    def update(self, problem, t, u, p):
        "Update problem at time t"

        # Add to accumulated CPU time
        timestep_cputime = time() - self._time
        self._cputime += timestep_cputime

        # Compute divergence
        if self.options["compute_divergence"]:
            check_divergence(u, p.function_space())

        # Update problem
        problem.update_problem(t, u, p)

        # Evaluate functional and error
        m = problem.reference(t)
        M = problem.functional(t, u, p)
        if m is None:
            e = None
            print "M = %g (missing reference value)" % M
        else:
            e = abs(M - m)
            print "M = %g (reference %g), error = %g (maximum %g)" % (M, m, e, max([e] + self._e))

        # Store values
        self._t.append(t)
        self._M.append(M)
        self._m.append(m)
        self._e.append(e)

        # save
        self.save(problem, t, u, p)


        # Check memory usage
        if self.options["check_mem_usage"]:
            if (self._timestep - 1) % self.options["check_frequency"] == 0:
                print 'Memory usage is:', self.getMyMemoryUsage()

        # Print progress
        print ""
        s = "Time step %d finished in %g seconds, %g%% done (t = %g, T = %g)." \
            % (self._timestep, timestep_cputime, 100.0 * (t / problem.T), t, problem.T)
        print s + "\n" + len(s) * "-"

        # Increase time step and record current time
        self._timestep += 1
        self._time = time()


    def eval(self):
        """"Return last functional value and maximum error in functional value on [0, T]"""

        # Return value
        if self._e[0] is None:
            return self._M[-1], None
        else:
            return self._M[-1], max([0.0] + self._e)

    def cputime(self):
        """"Return accumulated CPU time."""
        return self._cputime


def epsilon(u):
    """Return symmetric gradient."""
    return 0.5 * (grad(u) + grad(u).T)


def sigma(u, p, nu):
    """Return stress tensor."""
    return 2 * nu * epsilon(u) - p * Identity(u.cell().d)


def is_periodic(bcs):
    "Check if boundary conditions are periodic."
    # return all(isinstance(bc, PeriodicBC) for bc in bcs)
    return all(isinstance(bc, PeriodicBoundaryComputation) for bc in bcs)


def check_divergence(u, Q):
    """Check divergence of velocity."""

    # Compute L2 norm of divergence
    print "||div u||_L2 =", norm(u, "Hdiv0")

    # Compute projection of div u into Q_0
    pdivu = project(div(u), Q)
    zero = Constant(0.0)
    bc = DirichletBC(Q, zero, DomainBoundary())
    bc.apply(pdivu.vector())

    # Compute "weak" L2 norm of divergence
    print "||div u||_w  =", sqrt(abs(assemble(pdivu * div(u) * dx, mesh=Q.mesh())))
