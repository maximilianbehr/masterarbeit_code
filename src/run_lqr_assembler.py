from src.outputhandler import KarmanOutputHandler
from src.outputhandler import ProblemSolverOutputHandler
from src.outputhandler import LQRAssemblerOutputHandler
from src.lqr.lqr_assembler import LQR_Assembler


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
    "mat": None
}


# karman
REs = [100, 200, 300, 400, 500]
refinements = [2]
for RE in REs:
    for refinement in refinements:
        kohandler = KarmanOutputHandler()
        psohandler = ProblemSolverOutputHandler("karman", "stat_newton")
        lqrohandler = LQRAssemblerOutputHandler()
        OPTIONS["ref"] = refinement
        OPTIONS["RE"] = RE
        OPTIONS["mesh"] = kohandler.karman_mesh_xml(refinement)
        OPTIONS["boundaryfunction"] = kohandler.karman_boundary_xml(refinement)
        OPTIONS["u_stat"] = psohandler.u_xml(refinement, RE)
        OPTIONS["M_mtx"] = lqrohandler.M_mtx(refinement, RE)
        OPTIONS["S_mtx"] = lqrohandler.S_mtx(refinement, RE)
        OPTIONS["Mlower_mtx"] = lqrohandler.Mlower_mtx(refinement, RE)
        OPTIONS["Mupper_mtx"] = lqrohandler.Mupper_mtx(refinement, RE)
        OPTIONS["K_mtx"] = lqrohandler.K_mtx(refinement, RE)
        OPTIONS["R_mtx"] = lqrohandler.R_mtx(refinement, RE)
        OPTIONS["G_mtx"] = lqrohandler.G_mtx(refinement, RE)
        OPTIONS["Gt_mtx"] = lqrohandler.Gt_mtx(refinement, RE)
        OPTIONS["Blower_mtx"] = lqrohandler.Blower_mtx(refinement, RE)
        OPTIONS["Bupper_mtx"] = lqrohandler.Bupper_mtx(refinement, RE)
        OPTIONS["B_mtx"] = lqrohandler.B_mtx(refinement, RE)
        OPTIONS["C_mtx"] = lqrohandler.C_mtx(refinement, RE)
        OPTIONS["mat"] = lqrohandler.mat(refinement, RE)
        OPTIONS["options_json"] = lqrohandler.options_json(refinement, RE)

        lqrassembler = LQR_Assembler(OPTIONS)
        lqrassembler.unparameterized_lns_variational()
        lqrassembler.unparameterized_lns_ublas()
        lqrassembler.unparameterized_lns_npsc()
        lqrassembler.unparameterized_lns_mtx()
        lqrassembler.unparametrized_lns_mat()
        lqrassembler.save_options()