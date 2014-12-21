import os
import traceback

from outputhandler import KarmanOutputHandler
from outputhandler import ProblemSolverOutputHandler
from call_solver import call
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
           "T": None,
           "u_pvd": None,
           "p_pvd": None,
           "u_xml": None,
           "p_xml": None,
           "u_t_xml": None,
           "p_t_xml": None,
           "save_solution_at_t=T": True,
           "options_json": None,
           "log": None,
           "save_frequency": 10,
           "check_mem_usage": True,
           "check_frequency": 10,
           "compute_divergence": True,
           "debug": True,
           "dt": None,
           "dt_division": 0,
           "max_steps": None,
           "krylov_solver_absolute_tolerance": 1e-25,
           "krylov_solver_relative_tolerance": 1e-12,
           "krylov_solver_monitor_convergence": False,
           "form_compiler_optimize": True,
           "form_compiler_cpp_optimize": True
}


# karman
instant_clean = False
REs = [1, 2, 3, 4, 5, 10, 20, 50, 75, 100, 200, 300, 400, 500, 600, 700, 750]
refinements = [1, 2, 3, 4]

problem = "karman"

solvers = ["ipcs", "chorin", "instat_newton"]


for refinement in refinements:
    for solver in solvers:
        for RE in REs:
            kohandler = KarmanOutputHandler()
            psohandler = ProblemSolverOutputHandler(problem, solver, refinement, RE)

            OPTIONS["mesh"] = kohandler.karman_mesh_xml(refinement)
            OPTIONS["RE"] = RE
            OPTIONS["T"] = 12.0
            OPTIONS["u_pvd"] = psohandler.u_pvd()
            OPTIONS["p_pvd"] = psohandler.p_pvd()
            OPTIONS["u_xml"] = psohandler.u_xml()
            OPTIONS["p_xml"] = psohandler.p_xml()
            OPTIONS["u_t_xml"] = psohandler.u_t_xml()
            OPTIONS["p_t_xml"] = psohandler.p_t_xml()
            OPTIONS["options_json"] = psohandler.options_json()
            OPTIONS["log"] = psohandler.log()

            th = TeeHandler(OPTIONS["log"])
            th.start()

            try:
                if instant_clean:
                    os.system("instant-clean")

                print "{0:s}: Karman instat solver".format(gettime())
                call(problem, solver, OPTIONS)
                print "{0:s}: info Karman instat solver, ref={1:d}, RE={2:d}, solver={3:s} finished".format(
                    gettime(), refinement, RE, solver)
                th.stop()
            except Exception:
                traceback.print_exc()
                print "{0:s}: info Karman instat solver, ref={1:d}, RE={2:d}, solver={3:s} failed".format(gettime(),
                                                                                                          refinement,
                                                                                                          RE,
                                                                                                          solver)
                th.stop()
                break



