from pycmess import options
from pycmess import equation_dae2
from pycmess import PYCMESS_OP_TRANSPOSE
from pycmess import lrnm
from scipy.io import mmread
from scipy.io import mmwrite
from numpy import savetxt
import os
import json


class LQR_Solver():
    def __init__(self, argoptions):
        self.options = None
        self.opt = None
        self.eqn = None
        self.Z = None
        self.res2 = None

        self.options = argoptions.copy()

        # read data
        M = mmread(self.options["M_mtx"])
        S = mmread(self.options["S_mtx"])
        Mlower = mmread(self.options["Mlower_mtx"])
        Mupper = mmread(self.options["Mupper_mtx"])
        K = mmread(self.options["K_mtx"])
        R = mmread(self.options["R_mtx"])
        G = mmread(self.options["G_mtx"])
        Gt = mmread(self.options["Gt_mtx"])
        B = mmread(self.options["B_mtx"])
        C = mmread(self.options["C_mtx"])

        #setup equation
        self.eqn = equation_dae2()
        self.eqn.M = M
        self.eqn.A = -S - Mlower - Mupper - K - R
        self.eqn.G = -G
        self.eqn.B = B
        self.eqn.C = C
        self.eqn.delta = self.options["dae2_delta"]

    def setup_nm_adi_options(self):
        # setup nm and adi options
        self.opt = options()
        self.opt.adi.output = self.options["adi.output"]
        self.opt.nm.output = self.options["nm.output"]
        self.opt.nm.res2_tol = self.options["nm.res2_tol"]
        self.opt.adi.type = PYCMESS_OP_TRANSPOSE

    def solve(self):
        result = lrnm(self.eqn, self.opt)
        self.Z = result[0]
        self.res2 = result[1]
        # self.iter = result[2]


    def save(self):
        with open(self.options["Z_mtx"], "w") as handle:
            mmwrite(handle, self.Z)

        savetxt(self.options["res2_txt"], self.res2)

        fname = self.options["options_json"]
        if not os.path.exists(os.path.dirname(fname)):
            os.makedirs(os.path.dirname(fname))
        with open(fname, "w") as handle:
            json.dump(self.options, handle)

