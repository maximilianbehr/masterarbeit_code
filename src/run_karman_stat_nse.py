import os
import traceback
from call_solver import call
from outputhandler import KarmanOutputHandler
from outputhandler import ProblemSolverOutputHandler
from src.aux import gettime
from src.aux import TeeHandler
from src.aux import deletedir
from dolfin import parameters

# set dof reordering off
parameters["reorder_dofs_serial"] = False
parameters["form_compiler"]["optimize"] = True
parameters["form_compiler"].add("eliminate_zeros", True)

OPTIONS = {"mesh": None,
           "RE": None,
           "u_pvd": None,
           "p_pvd": None,
           "u_xml": None,
           "p_xml": None,
           "options_json": None,
           "logfile": None,
           "debug": False,
           "krylov_solver_absolute_tolerance": 1e-25,
           "krylov_solver_relative_tolerance": 1e-14,
           "krylov_solver_monitor_convergence": False,
           "form_compiler_optimize": True,
           "form_compiler_cpp_optimize": True,
           "newton_solver_max_iterations": 40,
           "newton_solver_absolute_tolerance": 1e-12,
           "newton_solver_relative_tolerance": 1e-14,
           "compute_divergence": True
}

# karman
instant_clean = False
REs = [1, 2, 3, 4, 5, 10, 20, 50, 75, 100, 200, 300, 400, 500, 600, 700, 750]
refinements = [1, 2, 3, 4]
# refinements = [5]


problem = "karman"
solver = "stat_newton"
for refinement in refinements:
    for RE in REs:

        kohandler = KarmanOutputHandler()
        psohandler = ProblemSolverOutputHandler(problem, solver, refinement, RE)

        OPTIONS["mesh"] = kohandler.karman_mesh_xml(refinement)
        OPTIONS["RE"] = RE
        OPTIONS["u_pvd"] = psohandler.u_pvd()
        OPTIONS["p_pvd"] = psohandler.p_pvd()
        OPTIONS["u_xml"] = psohandler.u_xml()
        OPTIONS["p_xml"] = psohandler.p_xml()
        OPTIONS["options_json"] = psohandler.options_json()
        OPTIONS["logfile"] = psohandler.log()

        th = TeeHandler(OPTIONS["logfile"])
        th.start()

        print "{0:s}: Karman Stationary ref = {1:d} RE = {2:d}".format(gettime(), refinement, RE)
        try:
            if instant_clean:
                os.system("instant-clean")

            call(problem, solver, OPTIONS)
            th.stop()
            print "--------------------------------------------------------------"
        except Exception, e:
            traceback.print_exc()
            th.stop()
            print "Solver Failed: Remove files and Directory"
            deletedir(psohandler.outputdir)
            # Reynoldsnumber to large increment refinement level
            break

