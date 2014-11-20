import os


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



#karman
instant_clean = False
REs = [100, 200, 300, 400, 500]
refinements = [3]
problems = ["karman"]
#solvers = ["chorin", "css1", "css2", "ipcs", "stat_newton"]
solvers = [ "stat_newton"]
for RE in REs:
    for refinement in refinements:
        for problem in problems:
            for solver in solvers:
                if instant_clean:
                    os.system("instant-clean")
                cmd = "python ns.py {0} {1} refinement_level={2} RE={3} debug=True".format(problem, solver, str(refinement),RE)
                print cmd
                os.system(cmd)
