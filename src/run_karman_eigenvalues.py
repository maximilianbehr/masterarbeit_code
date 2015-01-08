from src.outputhandler.lqrassembleroutputhandler import LQRAssemblerOutputHandler
from src.outputhandler.eigenoutputhandler import EigenOutputHandler
from src.aux import gettime
from src.aux import TeeHandler
from src.aux import deletedir
import os
from dolfin import parameters
from src.lqr.eigen import Eigen

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
    "Feed0_mtx": None,
    "dae2_delta": -0.02,
    "options_json": None,
    "logfile": None,
    "eig_eps": None,
    "eig_nopenalty_eps": None,
    "eig_bernoulli_eps": None,
    "eig_mtx": None,
    "eig_nopenalty_mtx": None,
    "eig_bernoulli_mtx": None
}

REs = [1, 2, 3, 4, 5, 10, 20, 50, 75, 100, 200]
refinements = [1]

for refinement in refinements:
    for RE in REs:

        lqrohandler = LQRAssemblerOutputHandler(refinement, RE)
        eighandler = EigenOutputHandler(refinement,RE)

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

        #take bernoulli feedback
        if os.path.isfile(lqrohandler.Feed0_mtx()):
            OPTIONS["Feed0_mtx"] = lqrohandler.Feed0_mtx()

        OPTIONS["options_json"] = eighandler.options_json()
        OPTIONS["logfile"] = eighandler.log()
        OPTIONS["eig_eps"] = eighandler.eig_eps()
        OPTIONS["eig_nopenalty_eps"] = eighandler.eig_nopenalty_eps()
        OPTIONS["eig_bernoulli_eps"] = eighandler.eig_bernoulli_eps()
        OPTIONS["eig_mtx"] = eighandler.eig_mtx()
        OPTIONS["eig_nopenalty_mtx"] = eighandler.eig_nopenalty_mtx()
        OPTIONS["eig_bernoulli_mtx"] = eighandler.eig_bernoulli_mtx()

        th = TeeHandler(OPTIONS["logfile"])
        th.start()

        try:

            eigensolver = Eigen(OPTIONS)

            print "{0:s}: Compute Eigenvalues".format(gettime())
            eigensolver.eigenvals()

            print "{0:s}: Compute Eigenvalues no penalty".format(gettime())
            eigensolver.eigenvals_nopenalty()

            if OPTIONS["Feed0_mtx"]:
                print "{0:s}: Compute Eigenvalues with Bernoulli Stabilization".format(gettime())
                eigensolver.eigenvals_bernoulli()
            th.stop()
        except Exception, e:
            print e
            th.stop()

