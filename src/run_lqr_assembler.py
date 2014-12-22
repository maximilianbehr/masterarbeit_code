from src.outputhandler import KarmanOutputHandler
from src.outputhandler import ProblemSolverOutputHandler
from src.outputhandler import LQRAssemblerOutputHandler
from src.lqr.lqr_assembler import LQR_Assembler
from src.aux import gettime
from src.aux import TeeHandler
from src.aux import deletedir
import traceback

from dolfin import parameters

# set dof reordering off
parameters["reorder_dofs_serial"] = False
parameters["form_compiler"]["optimize"] = True
parameters["form_compiler"].add("eliminate_zeros", True)


OPTIONS = {
    "ref": None,
    "RE": None,
    "mesh": None,
    "boundaryfunction": None,
    "u_stat": None,
    "penalty_eps": 0.001,
    "M_mtx": None,
    "S_mtx": None,
    "Mlower_mtx": None,
    "Mupper_mtx": None,
    "K_mtx": None,
    "R_mtx": None,
    "G_mtx": None,
    "Gt_mtx": None,
    "Blower_mtx": None,
    "Bupper_mtx": None,
    "B_mtx": None,
    "C_mtx": None,
    "mat": None,
    "options_json": None,
    "logfile": None,
    "observer_point1_x": 3.5,
    "observer_point1_y": 0.25,
    "observer_point2_x": 3.5,
    "observer_point2_y": 0.75
}


# karman
REs = [1, 2, 3, 4, 5, 10, 20, 50, 75, 100, 200]
refinements = [1, 2]
for refinement in refinements:
    for RE in REs:

        kohandler = KarmanOutputHandler()
        psohandler = ProblemSolverOutputHandler("karman", "stat_newton", refinement, RE)
        lqrohandler = LQRAssemblerOutputHandler(refinement, RE)
        OPTIONS["ref"] = refinement
        OPTIONS["RE"] = RE
        OPTIONS["mesh"] = kohandler.karman_mesh_xml(refinement)
        OPTIONS["boundaryfunction"] = kohandler.karman_boundary_xml(refinement)
        OPTIONS["u_stat"] = psohandler.u_xml()
        OPTIONS["M_mtx"] = lqrohandler.M_mtx()
        OPTIONS["S_mtx"] = lqrohandler.S_mtx()
        OPTIONS["Mlower_mtx"] = lqrohandler.Mlower_mtx()
        OPTIONS["Mupper_mtx"] = lqrohandler.Mupper_mtx()
        OPTIONS["K_mtx"] = lqrohandler.K_mtx()
        OPTIONS["R_mtx"] = lqrohandler.R_mtx()
        OPTIONS["G_mtx"] = lqrohandler.G_mtx()
        OPTIONS["Gt_mtx"] = lqrohandler.Gt_mtx()
        OPTIONS["Blower_mtx"] = lqrohandler.Blower_mtx()
        OPTIONS["Bupper_mtx"] = lqrohandler.Bupper_mtx()
        OPTIONS["B_mtx"] = lqrohandler.B_mtx()
        OPTIONS["C_mtx"] = lqrohandler.C_mtx()
        OPTIONS["mat"] = lqrohandler.mat()
        OPTIONS["options_json"] = lqrohandler.options_json_assembler()
        OPTIONS["logfile"] = lqrohandler.log_assembler()

        th = TeeHandler(OPTIONS["logfile"])
        th.start()
        try:
            print "{0:s}: LQR assembler refinement={1:d} RE={2:d}".format(gettime(), refinement, RE)
            lqrassembler = LQR_Assembler(OPTIONS)

            print "{0:s}: Build Variational Formulation".format(gettime())
            lqrassembler.lns_variational()

            print "{0:s}: Assemble uBLAS".format(gettime())
            lqrassembler.lns_ublas()

            print "{0:s}: Assemble NPSC".format(gettime())
            lqrassembler.lns_npsc()

            print "{0:s}: Save as mtx".format(gettime())
            lqrassembler.save_lns_mtx()

            print "{0:s}: Save as mat".format(gettime())
            lqrassembler.save_lns_mat()

            print "{0:s}: Save options".format(gettime())
            lqrassembler.save_options()

            th.stop()

        except Exception, e:
            traceback.print_exc()
            th.stop()
            print "Assembler Failed: Remove files and Directory"
            deletedir(lqrohandler.outputdir)
            break


