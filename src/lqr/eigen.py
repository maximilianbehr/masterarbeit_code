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
import pprint


class Eigen():
    def __init__(self, argoptions):

        #store and print options
        self.options = argoptions.copy()
        print pprint.pprint(self.options)
        fname = self.options["options_json"]

        if not os.path.exists(os.path.dirname(fname)):
            os.makedirs(os.path.dirname(fname))
        with open(fname, "w") as handle:
            json.dump(self.options, handle)


        #read and set matrices
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
        self.A = - self.S - self.Mlower - self.Mupper - self.K - self.R
        self.delta = self.options["dae2_delta"]



    def eigenvals(self):
        #build block matrices
        np = self.G.shape[1]
        upperblockM = hstack([self.M, self.delta*self.G])
        lowerblockM = hstack([self.delta*self.G.T, csr_matrix((np, np))])
        M = vstack([upperblockM, lowerblockM])
        upperblockA = hstack([self.A, self.G])
        lowerblockA = hstack([self.G.T, csr_matrix((np, np))])
        A = vstack([upperblockA, lowerblockA])

        #compute and sort eigenvalues by absolute value
        eigs = eigvals(A.todense(), M.todense(), overwrite_a=True, check_finite=False)
        eigs = eigs[numpy.argsort(numpy.absolute(eigs))]

        #write eigenvalues to file
        mmwrite(self.options["eig_mtx"], numpy.matrix(eigs))

        #split set of eigenvalues in stable, unstable and zeros (hopefully no)
        stable_eigs = eigs[eigs.real<0]
        unstable_eigs = eigs[eigs.real>0]
        zero_eigs = eigs[eigs.real==0]

        #plot eigenvalues
        fig, ax = plt.subplots()
        ax.plot(stable_eigs.real, stable_eigs.imag, "rx")
        ax.plot(unstable_eigs.real,unstable_eigs.imag, "bx")
        ax.plot(zero_eigs.real,zero_eigs.imag,"gx")
        plt.axvline(x=1.0/self.delta, linewidth=1, color="g",ls="dashed")
        xlimit = numpy.max(numpy.ceil(numpy.absolute(eigs.real)))
        ylimit = numpy.max(numpy.ceil(numpy.absolute(eigs.imag)))
        plt.xlim((-xlimit, xlimit))
        plt.ylim((-ylimit, ylimit))
        plt.xscale("symlog")
        plt.xlabel("Real")
        plt.ylabel("Imaginary")
        #plt.show()
        plt.savefig(self.options["eig_eps"])
        plt.close("all")

    def eigenvals_nopenalty(self):
        #build block matrices
        np = self.G.shape[1]
        upperblockM = hstack([self.M, self.delta*self.G])
        lowerblockM = hstack([self.delta*self.G.T, csr_matrix((np, np))])
        M = vstack([upperblockM, lowerblockM])
        upperblockA = hstack([-self.S - self.K - self.R, self.G])
        lowerblockA = hstack([self.G.T, csr_matrix((np, np))])
        A = vstack([upperblockA, lowerblockA])

        #compute and sort eigenvalues by absolute value
        eigs = eigvals(A.todense(), M.todense(), overwrite_a=True, check_finite=False)
        eigs = eigs[numpy.argsort(numpy.absolute(eigs))]

        #write eigenvalues to file
        mmwrite(self.options["eig_nopenalty_mtx"], numpy.matrix(eigs))

        #split set of eigenvalues in stable, unstable and zeros (hopefully no)
        stable_eigs = eigs[eigs.real<0]
        unstable_eigs = eigs[eigs.real>0]
        zero_eigs = eigs[eigs.real==0]

        #plot eigenvalues
        fig, ax = plt.subplots()
        ax.plot(stable_eigs.real, stable_eigs.imag, "rx")
        ax.plot(unstable_eigs.real,unstable_eigs.imag, "bx")
        ax.plot(zero_eigs.real,zero_eigs.imag,"gx")
        plt.axvline(x=1.0/self.delta, linewidth=1, color="g",ls="dashed")
        xlimit = numpy.max(numpy.ceil(numpy.absolute(eigs.real)))
        ylimit = numpy.max(numpy.ceil(numpy.absolute(eigs.imag)))
        plt.xlim((-xlimit, xlimit))
        plt.ylim((-ylimit, ylimit))
        plt.xscale("symlog")
        plt.xlabel("Real")
        plt.ylabel("Imaginary")
        #plt.show()
        plt.savefig(self.options["eig_nopenalty_eps"])
        plt.close("all")

    def eigenvals_bernoulli(self):

        if not self.options["Feed0_mtx"]:
            raise ValueError('No Bernoulli FeedBack given in options')
        else:
            self.Feed0 = mmread(self.options["Feed0_mtx"])

        #build block matrices
        np = self.G.shape[1]
        upperblockM = hstack([self.M, self.delta*self.G])
        lowerblockM = hstack([self.delta*self.G.T, csr_matrix((np, np))])
        M = vstack([upperblockM, lowerblockM])
        upperblockA = hstack([self.A-self.B*self.Feed0.T, self.G])
        lowerblockA = hstack([self.G.T, csr_matrix((np, np))])
        A = vstack([upperblockA, lowerblockA])

        #compute and sort eigenvalues by absolute value
        eigs = eigvals(A.todense(), M.todense(), overwrite_a=True, check_finite=False)
        eigs = eigs[numpy.argsort(numpy.absolute(eigs))]

        #write eigenvalues to file
        mmwrite(self.options["eig_bernoulli_mtx"], numpy.matrix(eigs))

        #split set of eigenvalues in stable, unstable and zeros (hopefully no)
        stable_eigs = eigs[eigs.real<0]
        unstable_eigs = eigs[eigs.real>0]
        zero_eigs = eigs[eigs.real==0]

        #plot eigenvalues
        fig, ax = plt.subplots()
        ax.plot(stable_eigs.real, stable_eigs.imag, "rx")
        ax.plot(unstable_eigs.real,unstable_eigs.imag, "bx")
        ax.plot(zero_eigs.real,zero_eigs.imag,"gx")
        plt.axvline(x=1.0/self.delta, linewidth=1, color="g",ls="dashed")
        xlimit = numpy.max(numpy.ceil(numpy.absolute(eigs.real)))
        ylimit = numpy.max(numpy.ceil(numpy.absolute(eigs.imag)))
        plt.xlim((-xlimit, xlimit))
        plt.ylim((-ylimit, ylimit))
        plt.xscale("symlog")
        plt.xlabel("Real")
        plt.ylabel("Imaginary")
        #plt.show()
        plt.savefig(self.options["eig_bernoulli_eps"])
        plt.close("all")
