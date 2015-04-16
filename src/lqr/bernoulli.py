import numpy as np
import scipy.linalg as scla
import scipy.io as scio
import scipy.sparse as scsp
import scipy.sparse.linalg as scspla
from dolfin import *
import os
from src.aux import createdir, write_matrix, gettime
import traceback


class Bernoulli():

    def __init__(self, const, ref, RE):

        # set parameters
        self.ref = ref
        self.RE = RE
        self.const = const
        self.bernoulli_instable_re = self.const.BERNOULLI_INSTABLE_RE


        # read compress system for simuation
        names = ["M", "M_BOUNDARY_CTRL", "S", "R", "K", "G", "GT", "B", "C"]
        self.mat = {}
        for name in names:
            self.mat[name] = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_NAME_MTX(ref, name, RE))
        if hasattr(self.mat["C"], 'toarray'):
            self.mat["C"] = self.mat["C"].toarray()

        # choose initial bernoulli strategy
        self.bernoullistrategy = self.const.BERNOULLI_STRATEGY_1

        # system sizes
        self.nv = self.mat["M"].shape[0]
        self.np = self.mat["G"].shape[1]

        # store blocks
        self.G = self.mat["G"]
        self.GT = self.mat["GT"]

        # build fullA
        upper = scsp.hstack([-self.mat["S"]-self.mat["R"]-self.mat["K"]-self.mat["M_BOUNDARY_CTRL"], self.mat["G"]])
        lower = scsp.hstack([self.mat["GT"], scsp.csr_matrix((self.np, self.np))])
        self.fullA = scsp.vstack([upper, lower])

        # build fullE
        upper = scsp.hstack([self.mat["M"], scsp.csr_matrix((self.nv, self.np))])
        lower = scsp.hstack([scsp.csr_matrix((self.np, self.nv)), scsp.csr_matrix((self.np, self.np))])
        self.fullE = scsp.vstack([upper, lower])

        # build B
        self.B = self.mat["B"]
        self.C = self.mat["C"]

        # set Feed0 a
        self.Feed0 = None
        self.Feed1 = None

    def _slepc_available(self):
        if not has_slepc():
            raise ValueError('FeNICS has no SLEPC support use scipy strategy instead')
        if not (has_petsc() and has_petsc4py()):
            raise ValueError('FeNICS has no PETSC support use scipy strategy instead')
        if not (has_linear_algebra_backend("PETSc")):
            raise ValueError("FeNICS has no PETSC linear algebra backend use scipy strategy instead")

    def _instable_subspace_shiftinvert(self, size):

        sigma = self.bernoullistrategy["sigma"]
        target = self.bernoullistrategy["target"]
        tol = self.bernoullistrategy["tol"]

        print "sigma={0:e} target={1:s} tol={2:e}".format(sigma, target, tol)
        if target == "LM" or target == "LR":
            # compute right eigenvectors
            print gettime(), "Compute right eigenvectors"
            try:
                self.dr, self.vr = scspla.eigs(self.fullA.tocsc(), size, self.fullE.tocsc(),  sigma=sigma, which=target, tol=tol)
            except scspla.ArpackNoConvergence as e:
                print "Arpack did not converge take current (right) eigenvalues eigenvectors"
                self.dr = e.eigenvalues
                self.vr = e.eigenvectors

            self.Ir = self.dr.real > 0
            self.nIr = self.Ir.sum()
            print "Instable Right Eigenvalues", self.dr[self.Ir]

            # compute left eigenvectors
            print gettime(), "Compute left eigenvectors"
            try:
                self.dl, self.vl = scspla.eigs(self.fullA.T.tocsc(), size, self.fullE.T.tocsc(),  sigma=sigma, which=target, tol=tol)
            except scspla.ArpackNoConvergence as e:
                print "Arpack did not converge take current (right) eigenvalues eigenvectors"
                self.dl = e.eigenvalues
                self.vl = e.eigenvectors

            self.Il = self.dl.real > 0
            self.nIl = self.Il.sum()
            print "Instable Left Eigenvalues", self.dl[self.Il]
        else:
            raise ValueError("target is propably not usefull for shiftinverting technique")

    def _instable_subspace_moebius(self, size):

        sigma = self.bernoullistrategy["sigma"]
        tau = self.bernoullistrategy["tau"]
        target = self.bernoullistrategy["target"]
        tol = self.bernoullistrategy["tol"]
        sortout = self.bernoullistrategy["sortout"]
        print "sigma={0:e}, tau={1:e}, tol={2:e}".format(sigma, tau, tol)

        # build linear operators for moebius transformed left and right ev problem
        shape = (self.nv+self.np, self.nv+self.np)
        Moeb1 = self.fullA-sigma*self.fullE
        Moeb2 = self.fullA-tau*self.fullE
        Moeb1 = Moeb1.tocsc(); Moeb2 = Moeb2.tocsc()
        Moeb1T = Moeb1.T.tocsc(); Moeb2T = Moeb2.T.tocsc()
        Moeb1solver = scspla.splu(Moeb1).solve; Moeb1Tsolver = scspla.splu(Moeb1T).solve
        matvecR = lambda x: Moeb1solver(Moeb2*x)
        matvecL = lambda x: Moeb1Tsolver(Moeb2T*x)
        MoebL = scspla.LinearOperator(shape, matvecL)
        MoebR = scspla.LinearOperator(shape, matvecR)

        print gettime(), "Compute right eigenvectors"
        try:
            dr, vr = scspla.eigs(MoebR, size, which=target, tol=tol)
        except scspla.ArpackNoConvergence as e:
            print "Arpack did not converge take current (right) eigenvalues eigenvectors"
            dr = e.eigenvalues
            vr = e.eigenvectors

        print gettime(), "Compute left eigenvectors"
        try:
            dl, vl = scspla.eigs(MoebL, size, which=target, tol=tol)
        except scspla.ArpackNoConvergence as e:
            print "Arpack did not converge take current (left) eigenvalues eigenvectors"
            dl = e.eigenvalues
            vl = e.eigenvectors

        # print found transformed eigenvalues
        print "Found (right/left) {0:d} {1:d} moebiues transformed eigenvalues ".format(dr.size, dl.size)
        print dr, dl

        # print transformed instable eigenvalues
        instabler = dr[np.absolute(dr)>1]
        instablel = dl[np.absolute(dl)>1]
        print "Found (right/left) {0:d} {1:d} moebiues transformed eigenvalues with |lambda|>1 (instable)".format(instabler.size, instablel.size)
        print instabler, np.absolute(instabler)
        print instablel, np.absolute(instablel)

        # transform instable eigenvalues back
        instablebackr = ((instabler*sigma-tau)/(instabler-1))
        instablebackl = ((instablel*sigma-tau)/(instablel-1))
        print "backtransformed (right/left) moebiues eigenvalues with |lambda|>1 (instable)"
        print instablebackr, instablebackl

        # sort out
        sortedout_r_idx = sortout(instablebackr)
        sortedout_l_idx = sortout(instablebackl)
        sortedout_r = instablebackr[sortedout_r_idx]
        sortedout_l = instablebackl[sortedout_l_idx]
        not_sortedout_r = instablebackr[np.logical_not(sortedout_r_idx)]
        not_sortedout_l = instablebackl[np.logical_not(sortedout_l_idx)]

        print "backtransformed  sorted out (right/left) moebiues eigenvalues with |lambda|>1 (instable)"
        print sortedout_r, sortedout_l

        print "backtransformed  not sorted out (right/left) moebiues eigenvalues with |lambda|>1 (instable)"
        print not_sortedout_r, not_sortedout_l

        # set attributes transform all back set attributes (eigenvectors change not)
        self.vr = vr
        self.vl = vl
        self.dr = ((dr*sigma-tau)/(dr-1))
        self.dl = ((dl*sigma-tau)/(dl-1))
        self.Ir = np.logical_and(self.dr.real > 0, np.logical_not(sortout(self.dr)))
        self.Il = np.logical_and(self.dl.real > 0, np.logical_not(sortout(self.dl)))
        #self.Ir = self.dr.real > 0
        #self.Il = self.dl.real > 0
        self.nIr = self.Ir.sum()
        self.nIl = self.Il.sum()

    def _resev(self):
        print "res right eigenvector"
        for i in np.nditer(np.where(self.Ir)):
            print "ev = {0:e}+{1:e}*i".format(self.dr[i].real, self.dr[i].imag)
            print "res = {0:e}".format(np.linalg.norm(self.fullA*self.vr[:, i]-self.dr[i]*self.fullE*self.vr[:, i]))

        print "res left eigenvector"
        for i in np.nditer(np.where(self.Il)):
            print "ev = {0:e}+{1:e}*i".format(self.dl[i].real, self.dl[i].imag)
            print "res = {0:e}".format(np.linalg.norm(self.fullA.T*self.vl[:, i]-self.dl[i]*self.fullE.T*self.vl[:, i]))

    def _instable_subspace_slepc(self, size):
        raise NotImplementedError('not further tested')
        # check if slepc is available
        self._slepc_available()

        # convert matrix to petsc
        from petsc4py import PETSc
        fullA = self.fullA.tocsr()
        fullE = self.fullE.tocsc()
        fullAT = self.fullA.T.tocsr()
        fullET = self.fullE.T.tocsc()

        fullApetsc = PETSc.Mat().createAIJ(size=fullA.shape, csr=(fullA.indptr, fullA.indices, fullA.data))
        fullEpetsc = PETSc.Mat().createAIJ(size=fullE.shape, csr=(fullE.indptr, fullE.indices, fullE.data))
        fullATpetsc = PETSc.Mat().createAIJ(size=fullAT.shape, csr=(fullAT.indptr, fullAT.indices, fullAT.data))
        fullETpetsc = PETSc.Mat().createAIJ(size=fullET.shape, csr=(fullET.indptr, fullET.indices, fullET.data))
        fullApetsc.assemble()
        fullEpetsc.assemble()
        fullATpetsc.assemble()
        fullETpetsc.assemble()

        fullApetscdolfin = PETScMatrix(fullApetsc)
        fullEpetscdolfin = PETScMatrix(fullEpetsc)
        fullATpetscdolfin = PETScMatrix(fullATpetsc)
        fullETpetscdolfin = PETScMatrix(fullETpetsc)

        eigensolver = SLEPcEigenSolver(fullApetscdolfin, fullEpetscdolfin)
        eigensolver.parameters["spectrum"] = "largest real"
        eigensolver.parameters["solver"] = "krylov-schur"
        eigensolver.parameters["problem_type"] = "pos_gen_non_hermitian"
        eigensolver.parameters["spectral_transform"] = "shift-and-invert"
        eigensolver.parameters["spectral_shift"] = self.const.BERNOULLI_STRATEGY["sigma"]
        eigensolver.parameters["verbose"] = self.const.BERNOULLI_STRATEGY["verbose"]
        eigensolver.parameters["maximum_iterations"] = self.const.BERNOULLI_STRATEGY["maxiter"]
        eigensolver.parameters["tolerance"] = self.const.BERNOULLI_STRATEGY["tol"]
        eps = eigensolver.eps()
        ksp = eps.getST().getKSP()
        ksp.setType("preonly")
        pc = ksp.getPC()
        pc.setType("lu")
        pc.setFactorSolverPackage("mumps")


        eigensolverT = SLEPcEigenSolver(fullATpetscdolfin, fullETpetscdolfin)
        eigensolverT.parameters["spectrum"] = "largest real"
        eigensolverT.parameters["solver"] = "krylov-schur"
        eigensolverT.parameters["problem_type"] = "pos_gen_non_hermitian"
        eigensolverT.parameters["spectral_transform"] = "shift-and-invert"
        eigensolverT.parameters["spectral_shift"] = self.const.BERNOULLI_STRATEGY["sigma"]
        eigensolverT.parameters["verbose"] = self.const.BERNOULLI_STRATEGY["verbose"]
        eigensolverT.parameters["maximum_iterations"] = self.const.BERNOULLI_STRATEGY["maxiter"]
        eigensolverT.parameters["tolerance"] = self.const.BERNOULLI_STRATEGY["tol"]
        epsT = eigensolverT.eps()
        kspT = epsT.getST().getKSP()
        kspT.setType("preonly")
        pcT = kspT.getPC()
        pcT.setType("lu")
        pcT.setFactorSolverPackage("mumps")

        # compute first tow desired eigenvalues
        eigensolver.solve(size)
        print "Number of Iterations     = {0:d}".format(eigensolver.get_iteration_number())
        print "Converged Eigenvalues    = {0:d}".format(eigensolver.get_number_converged())
        self.dr = np.zeros((eigensolver.get_number_converged(),), dtype=np.complex128)
        self.vr = np.zeros((self.np+self.nv, eigensolver.get_number_converged()), dtype=np.complex128)
        for i in range(0, eigensolver.get_number_converged()):
            r, c, rx, cx = eigensolver.get_eigenpair(i)
            ev = r+1j*c
            self.dr[i] = ev
            self.vr[:, i] = rx.array()+1j*cx.array()
            print "Right Eigenvalue {0:d} =".format(i), ev
        self.Ir = self.dr.real > 0
        self.nIr = self.Ir.sum()

        # compute first tow desired eigenvalues
        eigensolverT.solve(size)
        print "Number of Iterations     = {0:d}".format(eigensolverT.get_iteration_number())
        print "Converged Eigenvalues    = {0:d}".format(eigensolverT.get_number_converged())
        self.dl = np.zeros((eigensolverT.get_number_converged(),), dtype=np.complex128)
        self.vl = np.zeros((self.np+self.nv, eigensolverT.get_number_converged()), dtype=np.complex128)
        for i in range(0, eigensolverT.get_number_converged()):
            r, c, rx, cx = eigensolverT.get_eigenpair(i)
            ev = r+1j*c
            self.dl[i] = ev
            self.vl[:, i] = rx.array()+1j*cx.array()
            print "Left Eigenvalue {0:d} =".format(i), ev
        self.Il = self.dl.real > 0
        self.nIl = self.Il.sum()

    def _instable_subspace_fenics(self, size):
        raise NotImplementedError('This is no more supported')
        # check if slepc is available
        self._slepc_available()

        # build mesh and function spaces, stationary solution and dirichlet bcs and matrices
        nu = self.const.GET_NU(self.RE)
        penalty_eps = self.const.ASSEMBLER_PENALTY_EPS

        mesh = Mesh(self.const.MESH_XML(self.ref))
        V = VectorFunctionSpace(mesh, self.const.V, self.const.V_DIM)
        Q = FunctionSpace(mesh, self.const.Q, self.const.Q_DIM)
        W = V*Q
        boundaryfunction = MeshFunction("size_t", mesh, self.const.BOUNDARY_XML(self.ref))
        u_stat = Function(V, self.const.STATIONARY_U_XML(self.ref, self.RE))
        w_stat = Function(W, self.const.STATIONARY_W_XML(self.ref, self.RE))

        # get zero dirichlet boundary conditions
        noslip = Constant((0.0, 0.0))
        bcups = [DirichletBC(W.sub(0), noslip, boundaryfunction, INDICES) \
                      for INDICES in self.const.ASSEMBLER_COMPRESS_CTRL_DIRI_ZEROS]

        # build system matrices from weak formulation trial and test functions
        (dudt, dpdt) = TrialFunctions(W)
        (u, p) = TrialFunctions(W)
        (w_test, p_test) = TestFunctions(W)
        ds = Measure("ds")[boundaryfunction]
        M = inner(dudt, w_test) * dx; M_T = transpose(inner(dudt, w_test)) * dx
        S = nu * inner(grad(u), grad(w_test)) * dx; S_T = transpose(nu * inner(grad(u), grad(w_test))) * dx
        K = inner(grad(u_stat) * u, w_test) * dx; K_T = transpose(inner(grad(u_stat) * u, w_test)) * dx
        R = inner(grad(u) * u_stat, w_test) * dx; R_T = transpose(inner(grad(u) * u_stat, w_test) )*dx
        G = p * div(w_test) * dx; G_T = transpose(p * div(w_test)) * dx
        GT = div(u) * p_test * dx; GT_T = transpose(div(u) * p_test)*dx

        MGamma = sum([nu*(1.0/penalty_eps) * inner(u, w_test) * ds(indices)\
                for (controlfunc, indices, name) in self.const.ASSEMBLER_BOUNDARY_CONTROLS])

        MGamma_T = sum([transpose(nu*(1.0/penalty_eps) * inner(u, w_test))* ds(indices)\
                for (controlfunc, indices, name) in self.const.ASSEMBLER_BOUNDARY_CONTROLS])
        A = -1*(S+R+K)+G+GT -MGamma
        A_T = -1*(S_T+R_T+K_T)+G_T+GT_T -MGamma_T

        fullAdolfin = assemble(A)
        fullEdolfin = assemble(M)
        fullATdolfin = assemble(A_T)
        fullETdolfin = assemble(M_T)

        fullApetscdolfin = as_backend_type(fullAdolfin)
        fullEpetscdolfin = as_backend_type(fullEdolfin)
        fullATpetscdolfin = as_backend_type(fullATdolfin)
        fullETpetscdolfin = as_backend_type(fullETdolfin)

        # apply dirichlet bcs dont pertubate dirichlet nodes
        [bcup.apply(fullApetscdolfin) for bcup in bcups]
        [bcup.apply(fullEpetscdolfin) for bcup in bcups]
        [bcup.apply(fullATpetscdolfin) for bcup in bcups]
        [bcup.apply(fullETpetscdolfin) for bcup in bcups]
        PETScOptions().set("pc_type", "ilu")
        import ipdb
        ipdb.set_trace()
        eigensolver = SLEPcEigenSolver(fullApetscdolfin, fullEpetscdolfin)
        eigensolver.parameters["spectrum"] = self.const.BERNOULLI_STRATEGY["target"]
        eigensolver.parameters["solver"] = "krylov-schur"
        eigensolver.parameters["problem_type"] = "gen_non_hermitian"
        eigensolver.parameters["spectral_transform"] = self.const.BERNOULLI_STRATEGY["strategy"]
        eigensolver.parameters["spectral_shift"] = self.const.BERNOULLI_STRATEGY["sigma"]
        eigensolver.parameters["verbose"] = self.const.BERNOULLI_STRATEGY["verbose"]
        eigensolver.parameters["maximum_iterations"] = self.const.BERNOULLI_STRATEGY["maxiter"]
        #eigensolver.parameters["tolerance"] = self.const.BERNOULLI_STRATEGY["tol"]
        eigensolver.parameters["tolerance"] = 1e-5

        eigensolverT = SLEPcEigenSolver(fullATpetscdolfin, fullETpetscdolfin)
        eigensolverT.parameters["spectrum"] = self.const.BERNOULLI_STRATEGY["target"]
        eigensolverT.parameters["solver"] = "krylov-schur"
        eigensolverT.parameters["problem_type"] = "gen_non_hermitian"
        eigensolverT.parameters["spectral_transform"] = self.const.BERNOULLI_STRATEGY["strategy"]
        eigensolverT.parameters["spectral_shift"] = self.const.BERNOULLI_STRATEGY["sigma"]
        eigensolverT.parameters["verbose"] = self.const.BERNOULLI_STRATEGY["verbose"]
        eigensolverT.parameters["maximum_iterations"] = self.const.BERNOULLI_STRATEGY["maxiter"]
        eigensolverT.parameters["tolerance"] = self.const.BERNOULLI_STRATEGY["tol"]
        eigensolverT.parameters["tolerance"] = 1e-5

        # compute first tow desired eigenvalues
        eigensolver.solve(size)
        print "Number of Iterations     = {0:d}".format(eigensolver.get_iteration_number())
        print "Converged Eigenvalues    = {0:d}".format(eigensolver.get_number_converged())
        self.dr = np.zeros((eigensolver.get_number_converged(),), dtype=np.complex128)
        self.vr = np.zeros((self.np+self.nv, eigensolver.get_number_converged()), dtype=np.complex128)
        for i in range(0, eigensolver.get_number_converged()):
            r, c, rx, cx = eigensolver.get_eigenpair(i)
            ev = r+1j*c
            self.dr[i] = ev
            self.vr[:, i] = rx.array()+1j*cx.array()
            print "Right Eigenvalue {0:d} =".format(i), ev
        self.Ir = self.dr.real > 0
        self.nIr = self.Ir.sum()

        # compute first tow desired eigenvalues
        eigensolverT.solve(size)
        print "Number of Iterations     = {0:d}".format(eigensolverT.get_iteration_number())
        print "Converged Eigenvalues    = {0:d}".format(eigensolverT.get_number_converged())
        self.dl = np.zeros((eigensolverT.get_number_converged(),), dtype=np.complex128)
        self.vl = np.zeros((self.np+self.nv, eigensolverT.get_number_converged()), dtype=np.complex128)
        for i in range(0, eigensolverT.get_number_converged()):
            r, c, rx, cx = eigensolverT.get_eigenpair(i)
            ev = r+1j*c
            self.dl[i] = ev
            self.vl[:, i] = rx.array()+1j*cx.array()
            print "Left Eigenvalue {0:d} =".format(i), ev
        self.Il = self.dl.real > 0
        self.nIl = self.Il.sum()

    def _call_strategy(self):
        # chose strategy and solve
        size = self.bernoullistrategy["eigenvals"]
        strategy = self.bernoullistrategy["strategy"]
        solver = self.bernoullistrategy["solver"]

        if solver == "scipy" and strategy == "shiftinvert":
            self._instable_subspace_shiftinvert(size)
        elif solver == "scipy" and strategy == "moebius":
            self._instable_subspace_moebius(size)
        elif solver == "slepc":
            self._instable_subspace_slepc(size)
        elif solver == "fenics":
            self._instable_subspace_fenics(size)
        else:
            raise ValueError('unknown bernoulli eigenvalue shifting strategy {0:s}'.format(strategy))

    def solve(self):

        # call strategy
        self._call_strategy()

        if self.nIr != self.nIl or (self.bernoulli_instable_re >= self.RE and self.nIr==self.nIl):
            print "Number of left eigenvalues {0:d} != {1:d} number of right eigenvalues".format(self.nIl, self.nIr)
            print "Use Second Strategy"
            self.bernoullistrategy = self.const.BERNOULLI_STRATEGY_2
            self._call_strategy()
            if self.nIr != self.nIl or (self.bernoulli_instable_re >= self.RE and self.nIr==self.nIl):
                print "Use Third and last Strategy"
                self.bernoullistrategy = self.const.BERNOULLI_STRATEGY_3
                self._call_strategy()
                if self.nIr != self.nIl or (self.bernoulli_instable_re >= self.RE and self.nIr==self.nIl):
                    raise ValueError("Could not properly compute instable subspace.")

        elif self.nIl == 0 and self.nIr == 0:
            print "No instable eigenvalues detected, no need for initial feedback"
            return
        else:
            print "Found {0:d}, {1:d} instable eigenvalues".format(self.nIr, self.nIl)

        # print residuals
        self._resev()


        # sort instable eigenvektors and eigenvalues
        LEV = self.vl[:, self.Il]
        REV = self.vr[:, self.Ir]

        # project
        tA = np.dot(LEV.T, (self.fullA*REV))
        tE = np.dot(LEV.T, (self.fullE*REV))
        tB = np.dot(LEV.T, np.vstack((self.B, np.zeros((self.np, self.B.shape[1])))))
        tC = np.dot(np.hstack((self.C, np.zeros((self.C.shape[0], self.np)))), REV)

        # solve algebraic bernoulli equation on projected space
        try:
            XK, it = self._abe_gsign(tA, np.dot(tB, tB.T), tE)
            print "ABE solved within {0:d} steps.".format(it)

            # build initial feedback for nontransposed A
            K0 = LEV.T*self.fullE
            Blift = np.vstack((self.B, np.zeros((self.np, self.B.shape[1]))))
            K0 = np.dot(np.dot(np.dot(Blift.T, LEV), XK), K0)
            self.Feed0 = K0[:, 0:self.nv].real.T
            print "Computed Initial Feedback for transposed Riccatti equation"

        except Exception, e:
            traceback.print_exc(e)
            return

        # solve algebraic bernoulli equation on projected space
        try:
            XL, it = self._abe_gsign(tA.T, np.dot(tC.T, tC), tE.T)
            print "ABE solved within {0:d} steps.".format(it)

            # build initial feedback for transposed A
            K1 = REV.T*self.fullE.T
            Clift = np.hstack((self.C, np.zeros((self.C.shape[0], self.np))))
            K1 = np.dot(np.dot(np.dot(Clift, REV), XL), K1)
            self.Feed1 = K1[:, 0:self.nv].real.T
            print "Computed Initial Feedback for (non)transposed Riccatti equation"
        except Exception, e:
            traceback.print_exc(e)
            return


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
            write_matrix(filename, self.Feed0, "bernoulli Feed0, ref={0:d} RE={1:d}".format(self.ref, self.RE))

        if self.Feed1 is not None:
            filename = self.const.BERNOULLI_FEED1_CPS_MTX(self.ref, self.RE)
            createdir(filename)
            write_matrix(filename, self.Feed1, "bernoulli Feed1, ref={0:d} RE={1:d}".format(self.ref, self.RE))
