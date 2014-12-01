import os

from outputhandler import KarmanOutputHandler
from outputhandler import ProblemSolverOutputHandler
from call_solver import call


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
           "save_frequency": 10,
           "check_mem_usage": True,
           "check_frequency": 10,
           "compute_divergence": True,
           "debug": True,
           "dt": None,
           "krylov_solver_absolute_tolerance": 1e-25,
           "krylov_solver_relative_tolerance": 1e-12,
           "krylov_solver_monitor_convergence": False,
           "form_compiler_optimize": True,
           "form_compiler_cpp_optimize": True,
}




# karman
instant_clean = False
REs = [100,200,300,400,500,600,700]
refinements = [1,2,3,4,4,5]
problems = ["karman"]
solvers = ["chorin", "ipcs"]
for RE in REs:
    for refinement in refinements:
        for problem in problems:
            for solver in solvers:
                kohandler = KarmanOutputHandler()
                psohandler = ProblemSolverOutputHandler(problem, solver)

                OPTIONS["mesh"] = kohandler.karman_mesh_xml(refinement)
                OPTIONS["RE"] = RE
                OPTIONS["T"] = 20.0
                OPTIONS["u_pvd"] = psohandler.u_pvd(refinement, RE)
                OPTIONS["p_pvd"] = psohandler.p_pvd(refinement, RE)
                OPTIONS["u_xml"] = psohandler.u_xml(refinement, RE)
                OPTIONS["p_xml"] = psohandler.p_xml(refinement, RE)
                OPTIONS["u_t_xml"] = psohandler.u_t_xml(refinement, RE)
                OPTIONS["p_t_xml"] = psohandler.p_t_xml(refinement, RE)
                OPTIONS["options_json"] = psohandler.options_json(refinement, RE)

                if instant_clean:
                    os.system("instant-clean")
                call(problem, solver, OPTIONS)


