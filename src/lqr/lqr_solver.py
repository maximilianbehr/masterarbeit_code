from linearized_ctrl import LinearizedCtrl
import src.karman_const as const
from pycmess import equation_dae2, options, lrnm, PYCMESS_OP_TRANSPOSE
import numpy as np
import scipy.io as scio


class LQR_Solver():
    def __init__(self, ref, RE, RE_initial):

        # set parameters
        self.ref = ref
        self.RE = RE

        # get linearized matrices set arbitray values for dt, T, pertubationeps
        self.linearized = LinearizedCtrl(ref, RE, 1, 1, 1)

        # setup dae2 equation
        self.eqn = equation_dae2()
        self.eqn.M = self.linearized.Mcps
        self.eqn.A = - self.linearized.Scps \
                     - self.linearized.Kcps \
                     - self.linearized.Rcps \
                     - self.linearized.Mlowercps \
                     - self.linearized.Muppercps \

        self.eqn.G = self.linearized.Gcps
        self.eqn.B = self.linearized.Bcps
        self.eqn.C = self.linearized.Ccps
        self.eqn.delta = const.LINEARIZED_CTRL_DELTA

        # setup options
        self.opt = options()
        self.opt.nm.output = const.LINEARIZED_CTRL_NM_OUTPUT
        self.opt.nm.res2_tol = const.LINEARIZED_CTRL_NM_RES2
        self.opt.nm.rel2_change_tol = const.LINEARIZED_CTRL_NM_REL2_CHANGE
        self.opt.nm.rel_change_tol = const.LINEARIZED_CTRL_NM_REL_CHANGE
        self.opt.nm.maxit = const.LINEARIZED_CTRL_NM_MAXIT


        if RE_initial:
            self.opt.nm.nm_K0 = scio.mmread(const.LINEARIZED_CTRL_KINF_CPS_MTX(ref, RE_initial))

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

        #safe system
        self.linearized.save_compressed_matrices()

        #with open(const.LINEARIZED_CTRL_Z_CPS_MTX(self.ref, self.RE), "w") as handle:
        #    scio.mmwrite(handle, self.Z)

        with open(const.LINEARIZED_CTRL_KINF_CPS_MTX(self.ref, self.RE), "w") as handle:
            scio.mmwrite(handle, self.Kinf)
