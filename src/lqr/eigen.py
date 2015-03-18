import scipy.io as scio
import scipy.sparse as scsp
import scipy.sparse.linalg as scspla
import scipy.linalg as scla
from src.aux import *
import numpy as np
import matplotlib.pyplot as plt


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
            try:
                self.mat[name] = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_NAME_MTX(ref, name, RE))
            except:
                print "Could not load matrix {0:s}".format(name)
        # system sizes
        self.nv, self.np = self.mat["G"].shape
        print "nv={0:d} np={1:d}".format(self.nv, self.np)

        # build M
        upperblockM = scsp.hstack([self.mat["M"], const.LQR_DELTA*self.mat["G"]])
        lowerblockM = scsp.hstack([const.LQR_DELTA*self.mat["GT"], scsp.csr_matrix((self.np, self.np))])
        self.M = scsp.vstack([upperblockM, lowerblockM])

        # build A
        upperblockA = scsp.hstack([-self.mat["S"]-self.mat["R"]-self.mat["K"]-self.mat["M_BOUNDARY_CTRL"], self.mat["G"]])
        lowerblockA = scsp.hstack([self.mat["GT"], scsp.csr_matrix((self.np, self.np))])
        self.A = scsp.vstack([upperblockA, lowerblockA])

        # lift input matrix and feedback
        self.Bsys = np.vstack((self.mat["B"], np.zeros((self.np, self.mat["B"].shape[1]))))
        try:
            self.Kinfsys = np.vstack((self.mat["Kinf"], np.zeros((self.np, self.mat["Kinf"].shape[1]))))
        except:
            print "Probably no feedback available"

        # attributes for eigenvalues
        self.eig_sys = None
        self.eig_ric = None
        self.eig_ber = None

    @profile
    def compute_eig_sys(self):
        print gettime()
        # compute and sort eigenvalues by absolute value
        #if self.nv >= 10000:
        #A = self.A.tocsc()
        #M = self.M.tocsc()
        #self.eig_sys = scspla.eigs(A=A, k=400, M=M, which='LM', return_eigenvectors=False, sigma=0.0)
        #self.eig_sys = self.eig_sys[np.argsort(np.absolute(self.eig_sys))]

        #else:
        A = self.A.todense()
        M = self.M.todense()
        self.eig_sys = scla.eigvals(A, M, overwrite_a=True, check_finite=False)
        self.eig_sys = self.eig_sys[np.argsort(np.absolute(self.eig_sys))]

    @profile
    def compute_eig_ric(self):
        print gettime()
        # build A stable riccati feedback
        A_ric = self.A - np.dot(self.Bsys, self.Kinfsys.T)
        A_ric = A_ric.todense()

        # compute and sort eigenvalues by absolute value
        self.eig_ric = scla.eigvals(A_ric, self.M, overwrite_a=True, check_finite=False)
        self.eig_ric = self.eig_ric[np.argsort(np.absolute(self.eig_ric))]

    @profile
    def compute_eig_ber(self):
        print gettime()
        if not os.path.isfile(self.const.BERNOULLI_FEED0_CPS_MTX(self.ref, self.RE)):
            print "No initial Bernoulli Feedback found for ref = {0:d} and RE = {1:d}".format(self.ref, self.RE)
            return

        # build A bernoulli stable feedback
        Kbernoulli = scio.mmread(self.const.BERNOULLI_FEED0_CPS_MTX(self.ref, self.RE))
        Kbernoullisys = np.vstack((Kbernoulli, np.zeros((self.np, Kbernoulli.shape[1]))))
        A_ber = self.A - np.dot(self.Bsys, Kbernoullisys.T)
        A_ber = A_ber.todense()
        # compute and sort eigenvalues by absolute value
        self.eig_ber = scla.eigvals(A_ber, self.M, overwrite_a=True, check_finite=False)
        self.eig_ber = self.eig_ber[np.argsort(np.absolute(self.eig_ber))]

    def save(self):

        if self.eig_sys is not None:
            file = self.const.EIGEN_SYS_CPS_MTX(self.ref, self.RE)
            createdir(file)
            #scio.mmwrite(file, np.matrix(self.eig_sys))
            print self.eig_sys
            write_matrix(file, np.matrix(self.eig_sys), "Eigen Sys ref={0:d} RE={1:d}".format(self.ref, self.RE))

        if self.eig_ber is not None:
            file = self.const.EIGEN_BER_CPS_MTX(self.ref, self.RE)
            createdir(file)
            #scio.mmwrite(file, np.matrix(self.eig_ber))
            write_matrix(file, np.matrix(self.eig_ber), "Eigen Ber ref={0:d} RE={1:d}".format(self.ref, self.RE))

        if self.eig_ric is not None:
            file = self.const.EIGEN_RIC_CPS_MTX(self.ref, self.RE)
            createdir(file)
            #scio.mmwrite(file, np.matrix(self.eig_ric))
            write_matrix(file, np.matrix(self.eig_ric), "Eigen Ric ref={0:d} RE={1:d}".format(self.ref, self.RE))

    def _sort_eigs(self, eigs):
        stable_eigs = eigs[eigs.real < 0]
        unstable_eigs = eigs[eigs.real > 0]
        zero_eigs = eigs[eigs.real == 0]
        return stable_eigs, unstable_eigs, zero_eigs

    def _plot(self, eigs):
        # sort eigs
        stable_eigs, unstable_eigs, zero_eigs = self._sort_eigs(eigs)

        # plot eigenvalues

        fig, ax = plt.subplots()
        ax.plot(stable_eigs.real, stable_eigs.imag, "rx")
        ax.plot(unstable_eigs.real, unstable_eigs.imag, "bx")
        ax.plot(zero_eigs.real, zero_eigs.imag, "gx")
        xlimit = np.max(np.ceil(np.absolute(eigs.real)))
        ylimit = np.max(np.ceil(np.absolute(eigs.imag)))
        plt.xlim((-xlimit, xlimit))
        plt.ylim((-ylimit, ylimit))
        plt.xscale("symlog")
        plt.xlabel(u"Realteil")
        plt.ylabel(u"Imagin\u00e4rteil")
        plt.grid()
        # plt.show()


    def plot(self):
        olderr = np.geterr()
        np.seterr("ignore", "ignore", "ignore", "ignore", "ignore")

        if self.eig_ric is not None:
            self._plot(self.eig_ric)
            plt.savefig(self.const.EIGEN_RIC_CPS_PLOT(self.ref, self.RE, "png"))
            plt.savefig(self.const.EIGEN_RIC_CPS_PLOT(self.ref, self.RE, "eps"))
            plt.savefig(self.const.EIGEN_RIC_CPS_PLOT(self.ref, self.RE, "jpeg"))
            plt.close("all")

        if self.eig_ber is not None:
            self._plot(self.eig_ber)
            plt.savefig(self.const.EIGEN_BER_CPS_PLOT(self.ref, self.RE, "png"))
            plt.savefig(self.const.EIGEN_BER_CPS_PLOT(self.ref, self.RE, "eps"))
            plt.savefig(self.const.EIGEN_BER_CPS_PLOT(self.ref, self.RE, "jpeg"))
            plt.close("all")

        if self.eig_sys is not None:
            self._plot(self.eig_sys)
            plt.savefig(self.const.EIGEN_SYS_CPS_PLOT(self.ref, self.RE, "png"))
            plt.savefig(self.const.EIGEN_SYS_CPS_PLOT(self.ref, self.RE, "eps"))
            plt.savefig(self.const.EIGEN_SYS_CPS_PLOT(self.ref, self.RE, "jpeg"))
            plt.close("all")

        np.seterr(**olderr)
