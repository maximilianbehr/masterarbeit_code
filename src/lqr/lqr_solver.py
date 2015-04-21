import numpy as np
import scipy.io as scio
import os
from pycmess import equation_dae2, options, lrnm, PYCMESS_OP_TRANSPOSE
import warnings
from src.aux import createdir, write_matrix


class LQR_Solver():
    def __init__(self, const, ref, RE):

        # set parameters
        self.ref = ref
        self.RE = RE
        self.const = const

        # read compress system for simuation
        names = ["M", "M_BOUNDARY_CTRL", "S", "R", "K", "G", "GT", "B", "C"]
        self.mat = {}
        for name in names:
            self.mat[name] = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_NAME_MTX(ref, name, RE))


        # setup dae2 equation
        self.eqn = equation_dae2()
        self.eqn.M = self.mat["M"]
        self.eqn.A = -self.mat["S"] - self.mat["K"] - self.mat["R"] - self.mat["M_BOUNDARY_CTRL"]
        self.eqn.G = self.mat["G"]
        self.eqn.B = self.mat["B"]
        self.eqn.C = self.mat["C"]
        if hasattr(self.eqn.C, "todense"):
            self.eqn.C = self.eqn.C.todense()

        self.eqn.delta = const.LQR_DELTA

        # setup options
        self.opt = options()
        self.opt.nm.output = const.LQR_NM_OUTPUT
        self.opt.nm.res2_tol = const.LQR_NM_RES2
        self.opt.nm.rel_change_tol = const.LQR_NM_REL_CHANGE
        self.opt.nm.maxit = const.LQR_NM_MAXIT

        if os.path.isfile(const.BERNOULLI_FEED0_CPS_MTX(ref, RE)):
            print "take initial stabilizing feedback {0:s}".format(const.BERNOULLI_FEED0_CPS_MTX(ref, RE))
            self.opt.nm.nm_K0 = scio.mmread(const.BERNOULLI_FEED0_CPS_MTX(ref, RE))

        self.opt.adi.output = const.LQR_ADI_OUTPUT
        self.opt.adi.res2_tol = const.LQR_ADI_RES2
        self.opt.adi.maxit = const.LQR_ADI_MAXIT
        self.opt.adi.shifts.arp_m = const.LQR_ADI_ARP_M
        self.opt.adi.shifts.arp_p = const.LQR_ADI_ARP_P
        self.opt.adi.shifts.l0 = const.LQR_ADI_L0

        self.opt.adi.type = PYCMESS_OP_TRANSPOSE
        self.opt.adi.rel_change_tol = const.LQR_ADI_REL_CHANGE_TOL
        self.opt.adi.shifts.arp_m = const.LQR_ADI_ARP_M
        self.opt.adi.shifts.arp_p = const.LQR_ADI_ARP_P
        self.opt.adi.memory_usage = const.LQR_MEMORY_USAGE
        self.opt.adi.shifts.paratype = const.LQR_PARATYPE
        #self.opt.nm.gpStep = const.LQR_NM_GP
        #self.opt.adi.gpStep = const.LQR_GP
        # setup empty fields for solve
        self.Z = None
        self.res2 = None
        self.Kinf = None


    def solve(self):
        self.Z, self.status = lrnm(self.eqn, self.opt)
        self.res2 = self.status.res2_norms

        # self.iter = result[2]
        if self.res2[-1] > self.const.LQR_NM_RES2_SAVE:
            raise ValueError("No convergence")

        ZTM = self.eqn.M.T.dot(self.Z).T
        BTZ = np.dot(self.eqn.B.T, self.Z)
        self.Kinf = np.dot(BTZ, ZTM).T

    def save(self):

        #with open(const.LINEARIZED_CTRL_Z_CPS_MTX(self.ref, self.RE), "w") as handle:
        #    scio.mmwrite(handle, self.Z)

        file = self.const.ASSEMBLER_COMPRESS_CTRL_NAME_MTX(self.ref, "Kinf", self.RE)
        createdir(file)

        #with open(file, "w") as handle:
        #    scio.mmwrite(handle, self.Kinf)
        write_matrix(file, self.Kinf, "lqr solver Kinf, ref={0:d} RE={1:d}".format(self.ref, self.RE))

        try:
            with open(self.const.LQR_LOG(self.ref, self.RE),'w') as f:
                f.write(self.status.__str__())
                f.close()
        except:
            print "an error occured writing the logfile for lqr"