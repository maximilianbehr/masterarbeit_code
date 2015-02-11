import src.karman_const as const
import scipy.sparse as scsp
import scipy.io as scio
from dolfin import *
import numpy as np
import os
from src.aux import createdir


class Assembler():

    def __init__(self, ref, RE):

        # set parameters
        self.ref = ref
        self.RE = RE
        self.penalty_eps = const.ASSEMBLER_PENALTY_EPS

        # mesh and function spaces
        self.mesh = Mesh(const.MESH_XML(ref))
        self.V = VectorFunctionSpace(self.mesh, const.ASSEMBLER_V, const.ASSEMBLER_V_DIM)
        self.Q = FunctionSpace(self.mesh, const.ASSEMBLER_Q, const.ASSEMBLER_Q_DIM)
        self.boundaryfunction = MeshFunction("size_t", self.mesh, const.BOUNDARY_XML(ref))
        self.gupper = const.ASSEMBLER_UPPER_CONTROL
        self.glower = const.ASSEMBLER_LOWER_CONTROL
        self.u_stat = Function(self.V, const.STATIONARY_U_XML(ref, RE))

        # build variational and matrices
        self._lns_variational()
        self._lns_ublas()
        self._buildC()
        self._lns_npsc()

    def _get_sparsedata(self, m):
        """Function returns a scipy CSR Matrix for given Matrix (dolfin) M"""
        try:
            rowptr, colptr, vals = m.data()
            rows = m.size(0)
            cols = m.size(1)
            return scsp.csr_matrix((vals, colptr, rowptr), (rows, cols))
        except:
            raise ValueError(u"Implement _get_sparsedata for Matrix Datatype {0:s}".format(type(m)))

    def _lns_variational(self):
        # trial functions
        dudt = TrialFunction(self.V)
        u = TrialFunction(self.V)
        p = TrialFunction(self.Q)

        # test functions
        w_test = TestFunction(self.V)
        p_test = TestFunction(self.Q)

        # weak formulation of linearized navier stokes eqn
        ds = Measure("ds")[self.boundaryfunction]
        self.var = {}

        self.var["M"] = inner(dudt, w_test) * dx
        self.var["S"] = 1.0 / self.RE * inner(grad(u), grad(w_test)) * dx
        self.var["Bupper"] = 1.0 / self.penalty_eps * inner(self.gupper, w_test) * ds(
            const.GAMMA_BALL_CTRLUPPER_INDICES)

        self.var["Blower"] = 1.0 / self.penalty_eps * inner(self.glower, w_test) * ds(
            const.GAMMA_BALL_CTRLLOWER_INDICES)

        self.var["Mupper"] = 1.0 / self.penalty_eps * inner(u, w_test) * ds(
            const.GAMMA_BALL_CTRLUPPER_INDICES)

        self.var["Mlower"] = 1.0 / self.penalty_eps * inner(u, w_test) * ds(
            const.GAMMA_BALL_CTRLLOWER_INDICES)

        self.var["K"] = inner(grad(self.u_stat) * u, w_test) * dx
        self.var["R"] = inner(grad(u) * self.u_stat, w_test) * dx
        self.var["G"] = p * div(w_test) * dx
        self.var["Gt"] = div(u) * p_test * dx

    def _lns_ublas(self):
        if not self.var:
            raise ValueError("Call unparameterized_lns_variational to initialize attribute var.")

        # store backend
        la_backend = parameters["linear_algebra_backend"]

        # assemble matrices for PETSC backend
        parameters["linear_algebra_backend"] = "uBLAS"
        self.ublas = {}
        for key, val in self.var.iteritems():
            self.ublas[key] = assemble(val)

        # restore backend
        parameters["linear_algebra_backend"] = la_backend

    def _buildC(self):
        print "P1({0:.2f}, {1:.2f})\t P1({2:.2f}, {3:.2f})".format(const.ASSEMBLER_OBSERVER_POINT1_X, const.ASSEMBLER_OBSERVER_POINT1_Y, const.ASSEMBLER_OBSERVER_POINT2_X, const.ASSEMBLER_OBSERVER_POINT2_Y)
        p1 = Point(const.ASSEMBLER_OBSERVER_POINT1_X, const.ASSEMBLER_OBSERVER_POINT1_Y)
        p2 = Point(const.ASSEMBLER_OBSERVER_POINT2_X, const.ASSEMBLER_OBSERVER_POINT2_Y)
        points = [p1, p2]
        dists = {p1: float("inf"), p2: float("inf")}
        idxs = {p1: 0, p2: 0}
        #points = [p1]
        #dists = {p1: float("inf")}
        #idxs = {p1: 0}

        # find cell indices with midpoint has minimal distance
        for c in cells(self.mesh):
            for p in points:
                cdist = c.midpoint().distance(p)
                if cdist < dists[p]:
                    dists[p] = cdist
                    idxs[p] = c.index()

        # collect vertical dofs of cells with midpoint has minimal distance
        verticaldofs = np.empty(0, dtype=np.int32)
        for idx in idxs.values():
            verticaldofs = np.union1d(verticaldofs, self.V.sub(1).dofmap().cell_dofs(idx))

        #buid matrix
        C = scsp.lil_matrix((verticaldofs.size, self.V.dim()))
        j = 0
        for i in verticaldofs:
            C[j, i] = 1
            j += 1

        return C.todense()

    def _lns_npsc(self):
        if not self.ublas:
            raise ValueError("Call unparameterized_lns_ublas( to initialize attribute ublas.")

        self.npsc = {}
        for key, val in self.ublas.iteritems():
            if key == "Bupper" or key == "Blower":
                self.npsc[key] = val.array().reshape(val.array().size, 1)
            else:
                self.npsc[key] = self._get_sparsedata(val)

        self.npsc["B"] = np.concatenate((self.npsc["Bupper"], self.npsc["Blower"]), axis=1)

        self.npsc["C"] = self._buildC()

    def save_mtx(self):

        if not self.npsc:
            raise ValueError("Call unparameterized_lns_npsc to initialize attribute npsc.")

        createdir(os.path.dirname(const.ASSEMBLER_M_MTX(self.ref, self.RE)))

        with open(const.ASSEMBLER_M_MTX(self.ref, self.RE), "w") as handle:
            self.npsc["M"].eliminate_zeros()
            scio.mmwrite(handle, self.npsc["M"])

        with open(const.ASSEMBLER_S_MTX(self.ref, self.RE), "w") as handle:
            self.npsc["S"].eliminate_zeros()
            scio.mmwrite(handle, self.npsc["S"])

        with open(const.ASSEMBLER_MLOWER_MTX(self.ref, self.RE), "w") as handle:
            self.npsc["Mlower"].eliminate_zeros()
            scio.mmwrite(handle, self.npsc["Mlower"])

        with open(const.ASSEMBLER_MUPPER_MTX(self.ref, self.RE), "w") as handle:
            self.npsc["Mupper"].eliminate_zeros()
            scio.mmwrite(handle, self.npsc["Mupper"])

        with open(const.ASSEMBLER_K_MTX(self.ref, self.RE), "w") as handle:
            self.npsc["K"].eliminate_zeros()
            scio.mmwrite(handle, self.npsc["K"])

        with open(const.ASSEMBLER_R_MTX(self.ref, self.RE), "w") as handle:
            self.npsc["R"].eliminate_zeros()
            scio.mmwrite(handle, self.npsc["R"])

        with open(const.ASSEMBLER_G_MTX(self.ref, self.RE), "w") as handle:
            self.npsc["G"].eliminate_zeros()
            scio.mmwrite(handle, self.npsc["G"])

        with open(const.ASSEMBLER_GT_MTX(self.ref, self.RE), "w") as handle:
            self.npsc["Gt"].eliminate_zeros()
            scio.mmwrite(handle, self.npsc["Gt"])

        with open(const.ASSEMBLER_BLOWER_MTX(self.ref, self.RE), "w") as handle:
            scio.mmwrite(handle, self.npsc["Blower"])

        with open(const.ASSEMBLER_BUPPER_MTX(self.ref, self.RE), "w") as handle:
            scio.mmwrite(handle, self.npsc["Bupper"])

        with open(const.ASSEMBLER_B_MTX(self.ref, self.RE), "w") as handle:
            scio.mmwrite(handle, self.npsc["B"])

        with open(const.ASSEMBLER_C_MTX(self.ref, self.RE), "w") as handle:
            scio.mmwrite(handle, self.npsc["C"])

    def save_mat(self):

        if not self.npsc:
            raise ValueError("Call unparameterized_lns_npsc to initialize attribute npsc.")

        createdir(const.ASSEMBLER_MAT(self.ref, self.RE))

        with open(const.ASSEMBLER_MAT(self.ref, self.RE), "w") as handle:
            scio.savemat(handle, self.npsc, do_compression=True)













