import os

# Default options
OPTIONS = {"refinement_level": 0,
           "dt_division": 0,
           "save_solution": True,
           "save_frequency": 20,
           "check_mem_usage": True,
           "check_frequency": 10,
           "save_solution_at_t=T": True,
           "save_xml": True,
           "compute_stress": True,
           "compute_divergence": True,
           "debug": False,
           "max_steps": None,
           "krylov_solver_absolute_tolerance": 1e-25,
           "krylov_solver_relative_tolerance": 1e-12,
           "krylov_solver_monitor_convergence": False,
}


parameters["refinement_algorithm"] = "bisection"
if __name__=="__main__":

    #test problems from nsbench
    instant_clean = False
    refinements = [3]
    problems = ["beltrami", "drivencavity"]
    solvers = ["chorin", "css1", "css2", "ipcs"]

    for refinement in refinements:
        for problem in problems:
            for solver in solvers:
                if instant_clean:
                    os.system("instant-clean")
                print "python call_solver.py {0} {1} refinement_level={2} debug=True".format(problem, solver, str(refinement))
                os.system("python call_solver.py {0} {1} refinement_level={2} debug=True".format(problem, solver, str(refinement)))








