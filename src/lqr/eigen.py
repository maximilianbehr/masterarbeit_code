import scipy.io as scio
import scipy.sparse as scsp
import scipy.linalg as scla
from src.aux import *
import numpy as np
import matplotlib.pyplot as plt
from src.aux import createdir
import matplotlib


class Eigen():
    def __init__(self, const, ref, RE):

        # set parameters
        self.ref = ref
        self.RE = RE
        self.const = const

        # load compress system
        names = ["M", "M_BOUNDARY_CTRL", "S", "R", "K", "G", "GT", "B", "C", "Kinf"]
        self.mat = {}
        for name in names:
            self.mat[name] = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_NAME_MTX(ref, name, RE))

        # system sizes
        self.nv, self.np = self.mat["G"].shape

        # build M
        upperblockM = scsp.hstack([self.mat["M"], const.LQR_DELTA*self.mat["G"]])
        lowerblockM = scsp.hstack([const.LQR_DELTA*self.mat["GT"], scsp.csr_matrix((self.np, self.np))])
        self.M = scsp.vstack([upperblockM, lowerblockM]).todense()

        # build A
        upperblockA = scsp.hstack([-self.mat["S"]-self.mat["R"]-self.mat["K"]-self.mat["M_BOUNDARY_CTRL"], self.mat["G"]])
        lowerblockA = scsp.hstack([self.mat["GT"], scsp.csr_matrix((self.np, self.np))])
        self.A = scsp.vstack([upperblockA, lowerblockA]).todense()

        # lift input matrix and feedback
        self.Bsys = np.vstack((self.mat["B"], np.zeros((self.np, self.mat["B"].shape[1]))))
        self.Kinfsys = np.vstack((self.mat["Kinf"], np.zeros((self.np, self.mat["Kinf"].shape[1]))))

        # attributes for eigenvalues
        self.eig_sys = None
        self.eig_ric = None
        self.eig_ber = None

    @profile
    def compute_eig_sys(self):
        # compute and sort eigenvalues by absolute value
        self.eig_sys = scla.eigvals(self.A, self.M, overwrite_a=False, check_finite=False)
        self.eig_sys = self.eig_sys[np.argsort(np.absolute(self.eig_sys))]

    @profile
    def compute_eig_ric(self):
        # build A stable riccati feedback
        A_ric = self.A - np.dot(self.Bsys, self.Kinfsys.T)

        # compute and sort eigenvalues by absolute value
        self.eig_ric = scla.eigvals(A_ric, self.M, overwrite_a=True, check_finite=False)
        self.eig_ric = self.eig_ric[np.argsort(np.absolute(self.eig_ric))]

    @profile
    def compute_eig_ber(self):
        if not os.path.isfile(self.const.BERNOULLI_FEED0_CPS_MTX(self.ref, self.RE)):
            print "No initial Bernoulli Feedback found for ref = {0:d} and RE = {1:d}".format(self.ref, self.RE)
            return

        # build A bernoulli stable feedback
        Kbernoulli = scio.mmread(self.const.BERNOULLI_FEED0_CPS_MTX(self.ref, self.RE))
        Kbernoullisys = np.vstack((Kbernoulli, np.zeros((self.np,))))
        A_ber = self.A - np.dot(self.Bsys,Kbernoullisys.T)

        # compute and sort eigenvalues by absolute value
        self.eig_ber = scla.eigvals(A_ber, self.M, overwrite_a=True, check_finite=False)
        self.eig_ber = self.eig_ber[np.argsort(np.absolute(self.eig_ber))]

    def save(self):

        if self.eig_sys is not None:
            file = self.const.EIGEN_SYS_CPS_MTX(self.ref, self.RE)
            createdir(file)
            scio.mmwrite(file, np.matrix(self.eig_sys))

        if self.eig_ber is not None:
            file = self.const.EIGEN_BER_CPS_MTX(self.ref, self.RE)
            createdir(file)
            scio.mmwrite(file, np.matrix(self.eig_ber))

        if self.eig_ric is not None:
            file = self.const.EIGEN_RIC_CPS_MTX(self.ref, self.RE)
            createdir(file)
            scio.mmwrite(file, np.matrix(self.eig_ric))


