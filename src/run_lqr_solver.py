from src.outputhandler import LQRAssemblerOutputHandler
from src.lqr.lqr_solver import LQR_Solver
from src.aux import gettime
from src.aux import TeeHandler
from src.aux import deletedir
from dolfin import parameters

# set dof reordering off
parameters["reorder_dofs_serial"] = False

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
    "adi.output": 0,
    "nm.output": 1,
    "nm.res2_tol": 1e-9,
    "nm.rel_change_tol": 1e-14,
    "res2_txt": None,
    "options_json": None,
    "logfile": None
}

# REs = [1, 2, 3, 4, 5, 10, 20, 50, 75, 100, 200]
REs = [1, 2, 3, 4, 5]
#refinements = [1, 2, 3, 4, 5]
refinements = [1, 2, 3]

for refinement in refinements:
    for RE in REs:

        lqrohandler = LQRAssemblerOutputHandler(refinement, RE)
        OPTIONS["ref"] = refinement
        OPTIONS["RE"] = RE
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
        OPTIONS["res2_txt"] = lqrohandler.res2_txt()
        OPTIONS["logfile"] = lqrohandler.log_solver()

        th = TeeHandler(OPTIONS["logfile"])
        th.start()

        try:
            print "{0:s}: Setup pycmess equation".format(gettime())
            lqrsolver = LQR_Solver(OPTIONS)

            print "{0:s}: Setup pycmess options".format(gettime())
            lqrsolver.setup_nm_adi_options()

            print "{0:s}: Solve".format(gettime())
            lqrsolver.solve()

            print "{0:s}: Save Results".format(gettime())
            lqrsolver.save()

            th.stop()

        except Exception, e:
            print e
            print "Solver Failed"
            th.stop()
            # Reynoldsnumber to large increment refinement level
            break

