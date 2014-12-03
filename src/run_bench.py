import os

from outputhandler import KarmanOutputHandler, ProblemSolverOutputHandler
from call_solver import call
from dolfin import parameters


OPTIONS = {"refinement_level": None,
           "u_pvd": None,
           "p_pvd": None,
           "u_xml": None,
           "p_xml": None,
           "u_t_xml": None,
           "p_t_xml": None,
           "save_solution_at_t=T": True,
           "options_json": None,
           "save_frequency": 10,
           "check_mem_usage": True,
           "check_frequency": 10,
           "compute_divergence": True,
           "debug": True,
           "max_steps": None,
           "dt": None,
           "dt_division": 0,
           "krylov_solver_absolute_tolerance": 1e-25,
           "krylov_solver_relative_tolerance": 1e-12,
           "krylov_solver_monitor_convergence": False,
           "form_compiler_optimize": True,
           "form_compiler_cpp_optimize": True,
}




# karman
instant_clean = False
REs = [100]
refinements = [1,2,3]
problems = ["beltrami"]
solvers = ["chorin", "css1", "css2", "ipcs"]
for RE in REs:
    for refinement in refinements:
        for problem in problems:
            for solver in solvers:
                kohandler = KarmanOutputHandler()
                psohandler = ProblemSolverOutputHandler(problem, solver,refinement, RE)

                OPTIONS["refinement_level"] = refinement
                OPTIONS["u_pvd"] = psohandler.u_pvd()
                OPTIONS["p_pvd"] = psohandler.p_pvd()
                OPTIONS["u_xml"] = psohandler.u_xml()
                OPTIONS["p_xml"] = psohandler.p_xml()
                OPTIONS["u_t_xml"] = psohandler.u_t_xml()
                OPTIONS["p_t_xml"] = psohandler.p_t_xml()
                OPTIONS["options_json"] = psohandler.options_json()

                if instant_clean:
                    os.system("instant-clean")
                call(problem, solver, OPTIONS)


