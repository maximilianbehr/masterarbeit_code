# -*- coding: utf-8 -*-

import time
import pprint
from dolfin.cpp.common import set_log_active
from dolfin.cpp.common import list_timings
from dolfin import parameters
from problems import Problem
from solvers import Solver


def save_results(problem, solver, num_dofs, cputime, wct, functional, error):
    """Save results to file."""

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
    filename = "../results/results.log"
    file = open(filename, "a")
    file.write("{0:s}, {1:s}, {2:s}, {3:d}, {4:.15g}, {5:.15g}, {6:.15g},, {7:s}\n"
               .format(time.asctime(), problem, solver, num_dofs, cputime, wct, functional, str(error)))
    file.close()

def call(problem_name, solver_name, inoptions):

    #copy options
    options = inoptions.copy()


    # Set global DOLFIN parameters
    parameters["form_compiler"]["optimize"] = options["form_compiler_optimize"]
    parameters["form_compiler"]["cpp_optimize"] = options["form_compiler_cpp_optimize"]
    parameters["krylov_solver"]["absolute_tolerance"] = options["krylov_solver_absolute_tolerance"]
    parameters["krylov_solver"]["relative_tolerance"] = options["krylov_solver_relative_tolerance"]
    parameters["krylov_solver"]["monitor_convergence"] = options["krylov_solver_monitor_convergence"]

    # Set debug level
    set_log_active(options["debug"])

    # Create problem and solver
    problem = Problem(problem_name, options)
    solver = Solver(solver_name, options)
    print "Problem: " + str(problem)
    print "Solver:  " + str(solver)
    print pprint.pprint(options)



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
    save_results(problem, solver, num_dofs, solver.cputime(), wct, functional,  error)

