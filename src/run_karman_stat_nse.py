import os
from call_solver import call
from outputhandler import KarmanOutputHandler
from outputhandler import ProblemSolverOutputHandler


OPTIONS = {"mesh": None,
           "RE": None,
           "u_pvd": None,
           "p_pvd": None,
           "u_xml": None,
           "p_xml": None,
           "options_json": None,
           "debug": True,
           "krylov_solver_absolute_tolerance": 1e-25,
           "krylov_solver_relative_tolerance": 1e-14,
           "krylov_solver_monitor_convergence": True,
           "form_compiler_optimize": True,
           "form_compiler_cpp_optimize": True,
           "newton_solver_max_iterations": 40,
           "newton_solver_absolute_tolerance": 1e-14,
           "newton_solver_relative_tolerance": 1e-14
}




# karman
instant_clean = False
REs = [1, 5, 10, 20, 50, 75, 100, 200, 300, 400, 500]
refinements = [1, 2, 3, 4]
problem = "karman"
solver = "stat_newton"
for refinement in refinements:
    for RE in REs:
        print "Karman Stationary ref = {0:d} RE = {1:d}".format(refinement,RE)
        kohandler = KarmanOutputHandler()
        psohandler = ProblemSolverOutputHandler(problem, solver)

        OPTIONS["mesh"] = kohandler.karman_mesh_xml(refinement)
        OPTIONS["RE"] = RE
        OPTIONS["u_pvd"] = psohandler.u_pvd(refinement, RE)
        OPTIONS["p_pvd"] = psohandler.p_pvd(refinement, RE)
        OPTIONS["u_xml"] = psohandler.u_xml(refinement, RE)
        OPTIONS["p_xml"] = psohandler.p_xml(refinement, RE)
        OPTIONS["options_json"] = psohandler.options_json(refinement, RE)

        if instant_clean:
            os.system("instant-clean")

        try:
            call(problem, solver, OPTIONS)
        except:
            print "Solver Failed: Remove files and Directory"
            os.remove(OPTIONS["options_json"])
            os.rmdir(os.path.dirname(OPTIONS["options_json"]))

