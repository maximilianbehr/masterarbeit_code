from linearized_sim import Linearized

class Eigen():
    def __init__(self, ref, RE):

        #set parameters
        self.ref = ref
        self.RE = RE

        #generate compressed matrices via linearized
        self.linearized = Linearized(ref, RE, None, None, None)

        #get some parameters
        self.np

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
