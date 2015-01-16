import os
import traceback
from call_solver import call
from src.outputhandler.karmanoutputhandler import KarmanOutputHandler
from src.outputhandler.problemsolveroutputhandler import ProblemSolverOutputHandler
from src.aux import gettime
from src.aux import TeeHandler
from src.aux import deletedir
from dolfin import parameters

# set dof reordering off
parameters["reorder_dofs_serial"] = False
parameters["form_compiler"]["optimize"] = True

# problemname and solvername
problem = "karman"
solver = "stat_newton"

# setup options structure
OPTIONS = {"mesh": None,
           "RE": None,
           "u_pvd": None,
           "p_pvd": None,
           "u_xml": None,
           "p_xml": None,
           "options_json": None,
           "logfile": None,
           "debug": True,
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

# set Reynoldsnumbers and refinements
REs = [500, 1000, 1500, 1750, 2000, 2200, 2300, 2500]
refinements = [1, 2]

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

