import numpy as np
import scipy.linalg as scla
import scipy.io as scio
import scipy.sparse as scsp
import scipy.sparse.linalg as scspla
import os
from src.aux import createdir, write_matrix, gettime



class Bernoulli():

    def __init__(self, const, ref, RE):

        # set parameters
        self.ref = ref
        self.RE = RE
        self.const = const
        self.bereigenvalues = self.const.BERNOULLI_EIGENVALUES

        # read compress system for simuation
        names = ["M", "M_BOUNDARY_CTRL", "S", "R", "K", "G", "GT", "B", "C"]
        self.mat = {}
        for name in names:
            self.mat[name] = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_NAME_MTX(ref, name, RE))

        # system sizes
        self.nv = self.mat["M"].shape[0]
        self.np = self.mat["G"].shape[1]

        # build fullA
        upper = scsp.hstack([-self.mat["S"]-self.mat["R"]-self.mat["K"]-self.mat["M_BOUNDARY_CTRL"], self.mat["G"]])
        lower = scsp.hstack([self.mat["GT"], scsp.csr_matrix((self.np, self.np))])
        self.fullA = scsp.vstack([upper, lower])

        # build fullE
        upper = scsp.hstack([self.mat["M"], scsp.csr_matrix((self.nv, self.np))])
        self.fullE = scsp.vstack([upper, scsp.csr_matrix((self.np, self.nv+self.np))])

        # set Feed0 a
        self.Feed0 = None

    def _instable_subspace(self):

        size = self.bereigenvalues

        # compute right eigenvectors
        print gettime()
        print "Compute right eigenvectors"
        self.dr, self.vr = scspla.eigs(self.fullA.tocsc(), size, self.fullE.tocsc(),  sigma=0.25, which="LM")
        self.Ir = self.dr.real > 0
        self.nIr = self.Ir.sum()
        print "Instable Right Eigenvalues", self.dr[self.Ir]

        # compute left eigenvectors
        print gettime()
        print "Compute left eigenvectors"
        self.dl, self.vl = scspla.eigs(self.fullA.T.tocsc(), size, self.fullE.T.tocsc(), sigma=0.25, which="LM")
        self.Il = self.dl.real > 0
        self.nIl = self.Il.sum()
        print "Instable Left Eigenvalues", self.dl[self.Il]

    def solve(self):

        self._instable_subspace()

        while self.nIl != self.nIr:
            print "Number of left eigenvalues {0:d} != {1:d} number of right eigenvalues".format(self.nIl, self.nIr)
            self.bereigenvalues *= 2
            print "increase search space to {0:d}".format(self.bereigenvalues)

            if self.bereigenvalues >= self.fullA.shape[0]/2:
                raise ValueError("Could not properly compute instable subspace")
            else:
                self._instable_subspace()


        print "Found {0:d}, {1:d} instable eigenvalues".format(self.nIr, self.nIl)

        if self.nIl == 0 and self.nIr == 0:
            print "No instable eigenvalues detected, no need for initial feedback"
            return

        # sort instable eigenvektors and eigenvalues
        LEV = self.vl[:, self.Il]
        REV = self.vr[:, self.Ir]

        # project
        tA = np.dot(LEV.T, (self.fullA*REV))
        tE = np.dot(LEV.T, (self.fullE*REV))
        tB = np.dot(LEV.T, np.vstack((self.mat["B"], np.zeros((self.np, self.mat["B"].shape[1])))))

        # solve algebraic bernoulli equation on projected space
        try:
            XK, it = self._abe_gsign(tA, np.dot(tB, tB.T), tE)
            print "ABE solved within {0:d} steps.".format(it)
        except Exception, e:
            return

        # build initial feedback for nontransposed A
        K0 = LEV.T*self.fullE
        Blift = np.vstack((self.mat["B"], np.zeros((self.np, self.mat["B"].shape[1]))))
        K0 = np.dot(np.dot(np.dot(Blift.T, LEV), XK), K0)
        self.Feed0 = K0[:, 0:self.nv].real.T

    def _abe_gsign(self, A, B, E):

        if not A.shape == B.shape:
            raise ValueError('A and B must have the same number of rows and cols!')

        n = A.shape[0]
        # parameters for iteration
        it = 0
        maxit = self.const.BERNOULLI_MAXIT
        eps = np.finfo(np.float).eps
        tol = 10*n*np.sqrt(eps)
        Err = 1
        onemore = 0
        convergence = Err <= tol

        Enrm = np.linalg.norm(E, 1)
        PE, LE, UE = scla.lu(E)
        de = np.abs(np.diag(UE))
        if np.any(de < n*eps*Enrm):
            raise RuntimeWarning('E must not be singular.')
        de = np.power(de, 1.0/n)
        de = np.prod(de)

        while (it < maxit) and ((not convergence) or (convergence and (onemore < 2))):
            P, L, U = scla.lu(A)
            Ainv = np.dot(E, np.linalg.solve(U, np.linalg.solve(L, P)))
            # determinant scaling
            da = np.abs(np.diag(U))
            cs = np.float(np.prod(np.power(da, 1.0/n)))/de
            A = (A/cs + cs*np.dot(Ainv, E))/2
            B = (B/cs + cs*np.dot(np.dot(Ainv, B), Ainv.T))/2
            Err = np.linalg.norm(A-E, 1)/Enrm
            it += 1
            print "Step {0:d}, conv. crit. = {1:e}".format(it, Err)
            convergence = Err <= tol
            if convergence:
                onemore += 1

        X = np.linalg.lstsq(np.vstack((B, E.T-A.T)), np.vstack((np.linalg.solve(E.T, A.T).T+np.eye(n), np.zeros((n, n)))))[0]
        X = (X+X.T)/2

        if it == maxit and Err > tol:
            raise RuntimeError('ABE_SIGN: No convergence in {:d} iterations'.format(maxit))

        return X, it

    def save(self):
        if self.Feed0 is not None:
            filename = self.const.BERNOULLI_FEED0_CPS_MTX(self.ref, self.RE)
            createdir(filename)

            #with open(filename, "w") as handle:
            #    scio.mmwrite(handle, self.Feed0)
            write_matrix(filename, self.Feed0, "bernoulli Feed0, ref={0:d} RE={1:d}".format(self.ref, self.RE))
