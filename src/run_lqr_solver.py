from src.outputhandler import LQRAssemblerOutputHandler
from src.lqr.lqr_solver import LQR_Solver


OPTIONS = {
    "ref": None,
    "RE": None,
    "M_mtx": None,
    "S_mtx": None,
    "Mlower_mtx": None,
    "Mupper_mtx": None,
    "K_mtx": None,
    "R_mtx": None,
    "G_mtx": None,
    "Gt_mtx": None,
    "B_mtx": None,
    "C_mtx": None,
    "Z_mtx": None,
    "dae2_delta": -0.02,
    "adi.output": 1,
    "nm.output": 1,
    "nm.res2_tol": 1e-9,
    "res2_txt": None,
    "options_json":None
}


REs = [1, 5, 10, 20, 50, 75, 100, 200, 300, 400, 500]
refinements = [2]
for RE in REs:
    for refinement in refinements:
        print "LQR solver refinement = {0:d} RE = {1:d}".format(refinement, RE)
        lqrohandler = LQRAssemblerOutputHandler()
        OPTIONS["ref"] = refinement
        OPTIONS["RE"] = RE
        OPTIONS["M_mtx"] = lqrohandler.M_mtx(refinement, RE)
        OPTIONS["S_mtx"] = lqrohandler.S_mtx(refinement, RE)
        OPTIONS["Mlower_mtx"] = lqrohandler.Mlower_mtx(refinement, RE)
        OPTIONS["Mupper_mtx"] = lqrohandler.Mupper_mtx(refinement, RE)
        OPTIONS["K_mtx"] = lqrohandler.K_mtx(refinement, RE)
        OPTIONS["R_mtx"] = lqrohandler.R_mtx(refinement, RE)
        OPTIONS["G_mtx"] = lqrohandler.G_mtx(refinement, RE)
        OPTIONS["Gt_mtx"] = lqrohandler.Gt_mtx(refinement, RE)
        OPTIONS["B_mtx"] = lqrohandler.B_mtx(refinement, RE)
        OPTIONS["C_mtx"] = lqrohandler.C_mtx(refinement, RE)
        OPTIONS["Z_mtx"] = lqrohandler.Z_mtx(refinement, RE)
        OPTIONS["options_json"] = lqrohandler.options_json_solver(refinement, RE)
        OPTIONS["res2_txt"] = lqrohandler.res2_txt(refinement, RE)

        print "Setup pycmess equation"
        lqrsolver = LQR_Solver(OPTIONS)
        print "Setup pycmess options"
        lqrsolver.setup_nm_adi_options()
        print "Solve"
        lqrsolver.solve()
        print "Save Results"
        lqrsolver.save()
