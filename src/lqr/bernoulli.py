import numpy as np
import scipy as sc
import scipy.linalg as scla
import scipy.io as scio
import scipy.sparse as scsp
import scipy.sparse.linalg as scspla
import src.karman_const as const
import os

def abe_gsign(A, B, E):

    if not A.shape == B.shape:
        raise ValueError('A and B must have the same number of rows and cols!')

    n = A.shape[0]
    # parameters for iteration
    it = 0
    maxit = const.BERNOULLI_MAXIT
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

class Bernoulli():

    def __init__(self, ref, RE):

        # set parameters
        self.ref = ref
        self.RE = RE

        # read matrices
        M = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_M_MTX(self.ref, self.RE))
        S = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_S_MTX(self.ref, self.RE))
        R = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_R_MTX(self.ref, self.RE))
        K = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_K_MTX(self.ref, self.RE))
        Mlower = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_MLOWER_MTX(self.ref, self.RE))
        Mupper = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_MUPPER_MTX(self.ref, self.RE))
        G = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_G_MTX(self.ref, self.RE))
        GT = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_GT_MTX(self.ref, self.RE))
        self.B = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_B_MTX(self.ref, self.RE))
        self.C = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_C_MTX(self.ref, self.RE))

        self.nv, self.np = G.shape

        # build fullA
        upper = scsp.hstack([-S-R-K-Mlower-Mupper, G])
        lower = scsp.hstack([GT, scsp.csr_matrix((self.np, self.np))])
        self.fullA = scsp.vstack([upper, lower])

        # build fullE
        upper = scsp.hstack([M, scsp.csr_matrix((self.nv, self.np))])
        self.fullE = scsp.vstack([upper, scsp.csr_matrix((self.np, self.nv+self.np))])

        # set Feed0 and Feed1 attribute
        self.Feed0 = None
        # self.Feed1 = None

    def solve(self):

        # compute right eigenvectors
        print "Compute right eigenvectors"
        dr, vr = scspla.eigs(self.fullA.tocsc(), const.BERNOULLI_EIGENVALUES, self.fullE.tocsc(),  sigma=0)
        Ir = dr.real > 0
        nIr = Ir.sum()

        # compute left eigenvectors
        print "Compute left eigenvectors"
        dl, vl = scspla.eigs(self.fullA.T.tocsc(), const.BERNOULLI_EIGENVALUES, self.fullE.T.tocsc(),  sigma=0)
        Il = dl.real > 0
        nIl = Il.sum()

        print "Found {0:d}, {0:d} instable eigenvalues".format(nIr, nIl)

        if nIl == 0 and nIr == 0:
            print "No instable eigenvalues detected, no need for initial feedback"
            return

        # sort instable eigenvektors and eigenvalues
        LEV = vl[:, Il]
        REV = vr[:, Ir]

        # project
        tA = np.dot(LEV.T, (self.fullA*REV))
        tE = np.dot(LEV.T, (self.fullE*REV))
        tB = np.dot(LEV.T, np.vstack((self.B, np.zeros((self.np, self.B.shape[1])))))
        # tC = np.dot(np.hstack((self.C, np.zeros((self.C.shape[0], self.np)))), REV)

        # solve algebraic bernoulli equation on projected space
        try:
            XK, it = abe_gsign(tA, np.dot(tB, tB.T), tE)
            print "ABE solved within {0:d} steps.".format(it)
        except Exception, e:
            return

        # build initial feedback for nontransposed A
        K0 = LEV.T*self.fullE
        Blift = np.vstack((self.B, np.zeros((self.np, self.B.shape[1]))))
        K0 = np.dot(np.dot(np.dot(Blift.T, LEV), XK), K0)
        self.Feed0 = K0[:, 0:self.nv].real.T

        # solve algebraic bernoulli equation on projected space
        # try:
        #    XL, it = abe_gsign(tA.T, np.dot(tC.T, tC), tE.T)
        #    print "ABE solved within {0:d} steps.".format(it)
        # except Exception, e:
        #    self.Feed0 = None
        #    self.Feed1 = None
        #    return

        # build inital feedback for transposed A
        # L0 = REV.T*self.fullE.T
        # Clift = np.hstack((self.C, np.zeros((self.C.shape[0], self.np))))
        # L0 = np.dot(np.dot(np.dot(Clift, REV), XL), L0)

        # build feedback
        # self.Feed0 = K0[:, 0:self.nv].real.T
        # self.Feed1 = L0[:, 0:self.nv].real.T

    def save(self):
        if self.Feed0 is not None:
            if not os.path.exists(os.path.dirname(const.BERNOULLI_FEED0_CPS_MTX(self.ref, self.RE))):
                os.makedirs(os.path.dirname(const.BERNOULLI_FEED0_CPS_MTX(self.ref, self.RE)))

            with open(const.BERNOULLI_FEED0_CPS_MTX(self.ref, self.RE), "w") as handle:
                scio.mmwrite(handle, self.Feed0)
