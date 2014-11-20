# -*- coding: utf-8 -*-

import os
import time
from dolfin.cpp.common import set_log_active
from dolfin.cpp.common import list_timings
from dolfin import parameters
from problems import Problem, problems
from solvers import Solver, solvers

from outputhandler import ProblemSolverOutputHandler


def save_results(problem, solver, num_dofs, cputime, wct, functional, dt_division, error):
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
    file.write("%s, %s, %s, %d, %.15g, %.15g, %.15g, %s, %s\n" %
               (time.asctime(), problem, solver, num_dofs, cputime, wct, functional, str(dt_division), str(error)))
    file.close()

def ns(problem_name, solver_name, inoptions):


    options = inoptions.copy()


    # Set global DOLFIN parameters
    parameters["form_compiler"]["optimize"] = options["form_compiler_optimize"]
    parameters["form_compiler"]["cpp_optimize"] = options["form_compiler_cpp_optimize"]
    parameters["krylov_solver"]["absolute_tolerance"] = options["krylov_solver_absolute_tolerance"]
    parameters["krylov_solver"]["relative_tolerance"] = options["krylov_solver_relative_tolerance"]
    parameters["krylov_solver"]["monitor_convergence"] = options["krylov_solver_monitor_convergence"]


def main(args):
    """Parse command-line arguments and run solver"""


    #get outputdir
    psohandler = ProblemSolverOutputHandler(problem_name,solver_name)

    options["outputdir"]=os.path.join(psohandler.outputdir(),str(options.get("RE","")))




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

