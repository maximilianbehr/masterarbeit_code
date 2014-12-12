from src.outputhandler import LQRAssemblerOutputHandler
from src.outputhandler import KarmanOutputHandler
from src.outputhandler import ProblemSolverOutputHandler
from src.outputhandler import Linearized_NSE_SIM_OutputHandler

from src.lqr.linearized_nse import Linearized_NSE_SIM
from src.aux import gettime
from src.aux import TeeHandler
from src.aux import deletedir
from dolfin import parameters

# set dof reordering off
parameters["reorder_dofs_serial"] = False
parameters["form_compiler"]["optimize"] = True
parameters["form_compiler"].add("eliminate_zeros", True)

OPTIONS = {
    "ref": None,
    "RE": None,
    "mesh": None,
    "M_mtx": None,
    "S_mtx": None,
    "Mlower_mtx": None,
    "Mupper_mtx": None,
    "K_mtx": None,
    "R_mtx": None,
    "G_mtx": None,
    "Gt_mtx": None,
    "B_mtx": None,
    "options_json": None,
    "logfile": None,
    "dt": 0.1,
    "T": 12,
    "u_pvd": None,
    "u_t_xml": None,
    "save_frequency": 10,
    "pertubation_eps": 0.001,
    "u_stat": None

}

REs = [1]
refinements = [2]

for refinement in refinements:
    for RE in REs:
        kohandler = KarmanOutputHandler()
        lqrohandler = LQRAssemblerOutputHandler(refinement, RE)
        psohandler = ProblemSolverOutputHandler("karman", "stat_newton", refinement, RE)
        lnsesimhandler = Linearized_NSE_SIM_OutputHandler(refinement, RE)

        OPTIONS["ref"] = refinement
        OPTIONS["RE"] = RE
        OPTIONS["mesh"] = kohandler.karman_mesh_xml(refinement)
        OPTIONS["M_mtx"] = lqrohandler.M_mtx()
        OPTIONS["S_mtx"] = lqrohandler.S_mtx()
        OPTIONS["Mlower_mtx"] = lqrohandler.Mlower_mtx()
        OPTIONS["Mupper_mtx"] = lqrohandler.Mupper_mtx()
        OPTIONS["K_mtx"] = lqrohandler.K_mtx()
        OPTIONS["R_mtx"] = lqrohandler.R_mtx()
        OPTIONS["G_mtx"] = lqrohandler.G_mtx()
        OPTIONS["Gt_mtx"] = lqrohandler.Gt_mtx()
        OPTIONS["B_mtx"] = lqrohandler.B_mtx()
        OPTIONS["C_mtx"] = lqrohandler.C_mtx()
        OPTIONS["Z_mtx"] = lqrohandler.Z_mtx()
        OPTIONS["options_json"] = lqrohandler.options_json_solver()
        OPTIONS["logfile"] = lqrohandler.log_solver()
        OPTIONS["u_pvd"] = lnsesimhandler.u_pvd()
        OPTIONS["u_t_xml"] = lnsesimhandler.u_t_xml()


        th = TeeHandler(OPTIONS["logfile"])
        th.start()

        # try:
        print "{0:s}: Setup Linearized_NSE_SIM equation".format(gettime())
        linearized_nse_sim = Linearized_NSE_SIM(OPTIONS)

        print "{0:s}: Solve Linearized Navier Stokes".format(gettime())
        linearized_nse_sim.solve_ode()

        th.stop()

        # except Exception, e:
        #    print e
        #    print "Solver Failed"
        #    th.stop()
        #    # Reynoldsnumber to large increment refinement level
        #    break

