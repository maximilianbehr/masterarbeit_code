from pycmess import options
from pycmess import equation_dae2
from pycmess import PYCMESS_OP_TRANSPOSE
from pycmess import lrnm
from numpy import savetxt
import os
import json
from scipy.io import mmwrite
from scipy.io import mmread
from scipy.sparse import hstack
from scipy.sparse import vstack
from scipy.sparse import csr_matrix
from scipy.linalg import eigvals
import numpy
import matplotlib.pyplot as plt


class LQR_Solver():
    def __init__(self, argoptions):
        self.options = None
        self.opt = None
        self.eqn = None
        self.Z = None
        self.res2 = None

        self.options = argoptions.copy()

        # read data
        self.M = mmread(self.options["M_mtx"])
        self.S = mmread(self.options["S_mtx"])
        self.Mlower = mmread(self.options["Mlower_mtx"])
        self.Mupper = mmread(self.options["Mupper_mtx"])
        self.K = mmread(self.options["K_mtx"])
        self.R = mmread(self.options["R_mtx"])
        self.G = mmread(self.options["G_mtx"])
        self.Gt = mmread(self.options["Gt_mtx"])
        self.B = mmread(self.options["B_mtx"])
        self.C = mmread(self.options["C_mtx"])

        #setup equation
        self.eqn = equation_dae2()
        self.eqn.M = self.M
        self.eqn.A = -self.S - self.Mlower - self.Mupper - self.K - self.R
        self.eqn.G = self.G
        self.eqn.B = self.B
        self.eqn.C = self.C
        self.eqn.delta = self.options["dae2_delta"]

    def setup_nm_adi_options(self):
        # setup nm and adi options
        self.opt = options()
        self.opt.nm.output = self.options["nm.output"]
        self.opt.nm.res2_tol = self.options["nm.res2_tol"]
        self.opt.nm.rel_change_tol = self.options["nm.rel_change_tol"]
        self.opt.nm.maxit = self.options["nm.maxit"]

        self.opt.adi.output = self.options["adi.output"]
        self.opt.adi.res2_tol = self.options["adi.res2_tol"]
        self.opt.adi.maxit = self.options["adi.maxit"]
        self.opt.adi.type = PYCMESS_OP_TRANSPOSE

        self.opt.adi.shifts.arp_m = self.options["adi.shifts.arp_m"]
        self.opt.adi.shifts.arp_p = self.options["adi.shifts.arp_p"]
        self.opt.adi.shifts.l0 = self.options["adi.shifts.l0"]



    def solve(self):
        result = lrnm(self.eqn, self.opt)
        self.Z = result[0]
        self.res2 = result[1]
        # self.iter = result[2]
        if self.res2[-1] > self.opt.nm.res2_tol:
            raise ValueError("Newton ADI did not converge")


    def save(self):
        with open(self.options["Z_mtx"], "w") as handle:
            mmwrite(handle, self.Z)

        savetxt(self.options["res2_txt"], self.res2)

        fname = self.options["options_json"]
        if not os.path.exists(os.path.dirname(fname)):
            os.makedirs(os.path.dirname(fname))
        with open(fname, "w") as handle:
            json.dump(self.options, handle)

    def eigenvals(self):

        np = self.eqn.G.shape[1]
        upperblockM = hstack([self.eqn.M, self.eqn.delta*self.eqn.G])
        lowerblockM = hstack([self.eqn.delta*self.eqn.G.T, csr_matrix((np, np))])
        M = vstack([upperblockM, lowerblockM])
        upperblockA = hstack([self.eqn.A, self.eqn.G])
        lowerblockA = hstack([self.eqn.G.T, csr_matrix((np, np))])
        A = vstack([upperblockA,lowerblockA])
        eigs = eigvals(A.todense(), M.todense(), overwrite_a=True, check_finite=False)
        eigs = eigs[numpy.argsort(numpy.absolute(eigs))]
        mmwrite(self.options["eig_mtx"], numpy.matrix(eigs))

        stable_eigs = eigs[eigs.real<0]
        unstable_eigs = eigs[eigs.real>=0]

        fig, ax = plt.subplots()
        ax.plot(stable_eigs.real, stable_eigs.imag, "rx")
        ax.plot(unstable_eigs.real,unstable_eigs.imag, "bx")
        plt.axvline(x=1.0/self.eqn.delta, linewidth=1, color="g",ls="dashed")
        xlimit = numpy.max(numpy.ceil(numpy.absolute(eigs.real)))
        ylimit = numpy.max(numpy.ceil(numpy.absolute(eigs.imag)))
        plt.xlim((-xlimit, xlimit))
        plt.ylim((-ylimit, ylimit))
        plt.xscale("symlog")
        plt.xlabel("Real")
        plt.ylabel("Imaginary")
        #plt.show()
        plt.savefig(self.options["eig_eps"])

    def eigenvals_nopenalty(self):

        np = self.eqn.G.shape[1]
        upperblockM = hstack([self.eqn.M, self.eqn.delta*self.eqn.G])
        lowerblockM = hstack([self.eqn.delta*self.eqn.G.T, csr_matrix((np, np))])
        M = vstack([upperblockM, lowerblockM])
        upperblockA = hstack([-self.S - self.K - self.R, self.eqn.G])
        lowerblockA = hstack([self.eqn.G.T, csr_matrix((np, np))])
        A = vstack([upperblockA, lowerblockA])
        eigs = eigvals(A.todense(), M.todense(), overwrite_a=True, check_finite=False)
        eigs = eigs[numpy.argsort(numpy.absolute(eigs))]
        mmwrite(self.options["eig_nopenalty_mtx"], numpy.matrix(eigs))

        stable_eigs = eigs[eigs.real<0]
        unstable_eigs = eigs[eigs.real>=0]

        fig, ax = plt.subplots()
        ax.plot(stable_eigs.real, stable_eigs.imag, "rx")
        ax.plot(unstable_eigs.real,unstable_eigs.imag, "bx")
        plt.axvline(x=1.0/self.eqn.delta, linewidth=1, color="g",ls="dashed")
        xlimit = numpy.max(numpy.ceil(numpy.absolute(eigs.real)))
        ylimit = numpy.max(numpy.ceil(numpy.absolute(eigs.imag)))
        plt.xlim((-xlimit, xlimit))
        plt.ylim((-ylimit, ylimit))
        plt.xscale("symlog")
        plt.xlabel("Real")
        plt.ylabel("Imaginary")
        #plt.show()
        plt.savefig(self.options["eig_nopenalty_eps"])









