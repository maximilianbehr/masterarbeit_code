import scipy.io as scio
import scipy.sparse as scsp
import scipy.linalg as scla
import numpy as np
import os
from src.aux import createdir

class Eigen():
    def __init__(self, const, ref, RE):

        # set parameters
        self.ref = ref
        self.RE = RE
        self.const = const

        # load compress system
        Mcps = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_M_MTX(self.ref, self.RE))
        Mlowercps = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_MLOWER_MTX(self.ref, self.RE))
        Muppercps = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_MUPPER_MTX(self.ref, self.RE))
        Scps = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_S_MTX(self.ref, self.RE))
        Rcps = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_R_MTX(self.ref, self.RE))
        Kcps = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_K_MTX(self.ref, self.RE))
        Gcps = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_G_MTX(self.ref, self.RE))
        GTcps = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_GT_MTX(self.ref, self.RE))
        Bcps = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_B_MTX(self.ref, self.RE))

        # system sizes
        self.nv, self.np = Gcps.shape

        # build M
        upperblockM = scsp.hstack([Mcps, const.LINEARIZED_CTRL_DELTA*Gcps])
        lowerblockM = scsp.hstack([const.LINEARIZED_CTRL_DELTA*GTcps, scsp.csr_matrix((self.np, self.np))])
        self.M = scsp.vstack([upperblockM, lowerblockM]).todense()

        # build A
        upperblockA = scsp.hstack([-Scps-Rcps-Kcps-Muppercps-Mlowercps, Gcps])
        lowerblockA = scsp.hstack([GTcps, scsp.csr_matrix((self.np, self.np))])
        self.A = scsp.vstack([upperblockA, lowerblockA]).todense()

        # lift input matrix and feedback
        self.Bsys = np.vstack((Bcps, np.zeros((self.np, Bcps.shape[1]))))

        # attributes for eigenvalues
        self.eig_sys = None
        self.eig_ric = None
        self.eig_ber = None

    def compute_eig_sys(self):
        # compute and sort eigenvalues by absolute value
        self.eig_sys = scla.eigvals(self.A, self.M, overwrite_a=True, check_finite=False)
        self.eig_sys = self.eig_sys[np.argsort(np.absolute(self.eig_sys))]

    def compute_eig_ric(self):
        # build A stable riccati feedback
        Kinfcps = scio.mmread(self.const.LINEARIZED_CTRL_KINF_CPS_MTX(self.ref, self.RE))
        Kinfsys = np.vstack((Kinfcps, np.zeros((self.nump, Kinfcps.shape[1]))))
        A_ric = self.A - self.Bsys*Kinfsys.T

        # compute and sort eigenvalues by absolute value
        self.eig_ric = scla.eigvals(A_ric, self.M, overwrite_a=True, check_finite=False)
        self.eig_ric = self.eig_ric[np.argsort(np.absolute(self.eig_ric))]

    def compute_eig_ber(self):
        if not os.path.isfile(self.const.BERNOULLI_FEED0_CPS_MTX(self.ref, self.RE)):
            print "No initial Bernoulli Feedback found for ref = {0:d} and RE = {1:d}".format(self.ref, self.RE)
            return

        # build A bernoulli stable feedback
        Kbernoulli = scio.mmread(self.const.BERNOULLI_FEED0_CPS_MTX(self.ref, self.RE))
        Kbernoullisys = np.vstack((Kbernoulli, np.zeros((self.np,))))
        A_ber = self.A - self.Bsys*Kbernoullisys.T

        # compute and sort eigenvalues by absolute value
        self.eig_ber = scla.eigvals(A_ber.todense(), self.M, overwrite_a=True, check_finite=False)
        self.eig_ber = self.eig_ber[np.argsort(np.absolute(self.eig_ber))]
        self.eig_ber = self.eig_ber[np.argsort(np.absolute(self.eig_ber))]


    def save(self):

        if self.eig_sys:
            file = self.const.EIGEN_SYS_CPS_MTX(self.ref, self.RE)
            createdir(file)
            scio.mmwrite(file, np.matrix(self.eig_sys))

        if self.eig_ber:
            file = self.const.EIGEN_BER_CPS_MTX(self.ref, self.RE)
            createdir(file)
            scio.mmwrite(file, np.matrix(self.eig_ber))

        if self.eig_ric:
            file = self.const.EIGEN_RIC_CPS_MTX(self.ref, self.RE)
            createdir(file)
            scio.mmwrite(file, np.matrix(self.eig_ric))


