import os
import src.karman_const as const
from pycmess import equation_dae2, options, lrnm, PYCMESS_OP_TRANSPOSE
import numpy as np
import scipy.io as scio
import warnings

class LQR_Solver():
    def __init__(self, ref, RE):

        # set parameters
        self.ref = ref
        self.RE = RE

        # setup dae2 equation
        self.eqn = equation_dae2()
        self.eqn.M = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_M_MTX(self.ref, self.RE))
        self.eqn.A = - scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_S_MTX(self.ref, self.RE)) \
                     - scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_K_MTX(self.ref, self.RE)) \
                     - scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_R_MTX(self.ref, self.RE)) \
                     - scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_MLOWER_MTX(self.ref, self.RE)) \
                     - scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_MUPPER_MTX(self.ref, self.RE))

        self.eqn.G = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_G_MTX(self.ref, self.RE))
        self.eqn.B = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_B_MTX(self.ref, self.RE))
        self.eqn.C = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_C_MTX(self.ref, self.RE))
        self.eqn.delta = const.LINEARIZED_CTRL_DELTA

        # setup options
        self.opt = options()
        self.opt.nm.output = const.LINEARIZED_CTRL_NM_OUTPUT
        self.opt.nm.res2_tol = const.LINEARIZED_CTRL_NM_RES2
        self.opt.nm.rel2_change_tol = const.LINEARIZED_CTRL_NM_REL2_CHANGE
        self.opt.nm.rel_change_tol = const.LINEARIZED_CTRL_NM_REL_CHANGE
        self.opt.nm.maxit = const.LINEARIZED_CTRL_NM_MAXIT


        if os.path.isfile(const.BERNOULLI_FEED0_CPS_MTX(ref, RE)):
            self.opt.nm.nm_K0 = scio.mmread(const.BERNOULLI_FEED0_CPS_MTX(ref, RE))
        elif RE >= 400:
            warnings.warn("large reynoldsnumber and no initial feedback")

        self.opt.adi.output = const.LINEARIZED_CTRL_ADI_OUTPUT
        self.opt.adi.res2_tol = const.LINEARIZED_CTRL_ADI_RES2
        self.opt.adi.maxit = const.LINEARIZED_CTRL_ADI_MAXIT
        self.opt.adi.type = PYCMESS_OP_TRANSPOSE

        # setup empty fields for solve
        self.Z = None
        self.res2 = None
        self.Kinf = None


    def solve(self):
        result = lrnm(self.eqn, self.opt)
        self.Z = result[0]
        self.res2 = result[1]
        # self.iter = result[2]
        if self.res2[-1] > const.LINEARIZED_CTRL_NM_RES2_SAVE:
            raise ValueError("No convergence")

        ZTM = self.eqn.M.T.dot(self.Z).T
        BTZ = np.dot(self.eqn.B.T, self.Z)
        self.Kinf = np.dot(BTZ, ZTM).T

    def save(self):

        #with open(const.LINEARIZED_CTRL_Z_CPS_MTX(self.ref, self.RE), "w") as handle:
        #    scio.mmwrite(handle, self.Z)

        file = const.LINEARIZED_CTRL_KINF_CPS_MTX(self.ref, self.RE)
        if not os.path.exists(file):
            os.makedirs(os.path.dirname(file))

        with open(file, "w") as handle:
            scio.mmwrite(handle, self.Kinf)
