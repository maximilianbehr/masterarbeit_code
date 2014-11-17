# -*- coding: utf-8 -*-

import sys

from src import problems
import time
from dolfin import set_log_active, parameters, list_timings
from src.problems import Problem
from solvers import Solver, solvers






# Default options
OPTIONS = {"refinement_level": 0,
           "dt_division": 0,
           "save_solution": True,
           "save_frequency": 20,
           "check_mem_usage": True,
           "check_frequency": 10,
           "save_solution_at_t=T": True,
           "save_xml": True,
           "save_TimeSeries_u": False,
           "save_TimeSeries_p": False,
           "plot_functional": False,
           "compute_stress": True,
           "compute_divergence": True,
           "debug": False,
           "max_steps": None,
           "krylov_solver_absolute_tolerance": 1e-25,
           "krylov_solver_relative_tolerance": 1e-12,
           "krylov_solver_monitor_convergence": False,
           "newton_solver_max_iterations": 25
}
parameters["refinement_algorithm"] = "bisection"


def save_results(problem, solver, num_dofs, cputime, wct, functional, dt_division, error):
    "Save results to file."

    # Print summary
    print ""
    print "Problem    |", problem
    print "Solver     |", solver
    print "Unknowns   |", num_dofs
    print "CPU time   |", cputime
    print "WCT time   |", wct
    print "Overhead   |", wct - cputime
    print "Functional |", functional
    print "Error      |", error

    # Print DOLFIN summary
    set_log_active(True)
    list_timings()

    # Append to file
    filename = "results/results.log"
    file = open(filename, "a")
    file.write("%s, %s, %s, %d, %.15g, %.15g, %.15g, %s, %s\n" %
               (time.asctime(), problem, solver, num_dofs, cputime, wct, functional, str(dt_division), str(error)))
    file.close()


def usage():
    "Print usage"
    print """\
Usage: ns problem solver

Available problems:

%s

Available solvers:

%s
""" % ("\n".join("  " + p for p in problems),
       "\n".join("  " + s for s in solvers))


def main(args):
    "Parse command-line arguments and run solver"

    # Check arguments
    if not len(args) >= 2:
        usage()
        return 2

    # Get problem and solver
    problem_name, solver_name = args[:2]

    # Get options
    options = OPTIONS.copy()
    for arg in args[2:]:
        try:
            key, value = arg.split("=")
            try:
                options[key] = eval(value)
            except:
                options[key] = str(value)
        except:
            print "Warning: Unhandled command-line argument", arg

    # Set global DOLFIN parameters
    parameters["form_compiler"]["optimize"] = True
    parameters["form_compiler"]["cpp_optimize"] = True
    parameters["krylov_solver"]["absolute_tolerance"] = options["krylov_solver_absolute_tolerance"]
    parameters["krylov_solver"]["relative_tolerance"] = options["krylov_solver_relative_tolerance"]
    parameters["krylov_solver"]["monitor_convergence"] = options["krylov_solver_monitor_convergence"]

    # Set debug level
    set_log_active(options["debug"])

    # Set refinement level
    dt_division = str(options["dt_division"])

    # Create problem and solver
    problem = Problem(problem_name, options)
    solver = Solver(solver_name, options)
    print "Problem: " + str(problem)
    print "Solver:  " + str(solver)

    # Solve problem with solver
    wct = time.time()
    u, p = solver.solve(problem)

    # Compute elapsed time
    wct = time.time() - wct

    # Compute number of degrees of freedom
    num_dofs = u.vector().size() + p.vector().size()

    # Get functional value and error
    functional, error = solver.eval()

    # Save results
    save_results(problem, solver, num_dofs, solver.cputime(), wct, functional, dt_division, error)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
