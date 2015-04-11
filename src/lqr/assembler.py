import scipy.sparse as scsp
import scipy.io as scio
from dolfin import *
import numpy as np
import os
from src.aux import createdir, write_matrix


class Assembler():

    def __init__(self, const, ref, RE):

        # set parameters
        self.ref = ref
        self.RE = RE
        self.nu = const.GET_NU(RE)
        self.penalty_eps = const.ASSEMBLER_PENALTY_EPS
        self.const = const

        # mesh and function spaces
        self.mesh = Mesh(const.MESH_XML(ref))
        self.V = VectorFunctionSpace(self.mesh, const.V, const.V_DIM)
        self.Q = FunctionSpace(self.mesh, const.Q, const.Q_DIM)
        self.boundaryfunction = MeshFunction("size_t", self.mesh, const.BOUNDARY_XML(ref))
        self.u_stat = Function(self.V, const.STATIONARY_U_XML(ref, RE))

        # build variational and matrices
        self._lns_variational()
        self._lns_ublas()
        self._lns_npsc()


    def _get_sparsedata(self, m):
        """Function returns a scipy CSR Matrix for given Matrix (dolfin) M"""
        try:
            rowptr, colptr, vals = m.data()
            rows = m.size(0)
            cols = m.size(1)
            return scsp.csr_matrix((vals, colptr, rowptr), (rows, cols))
        except:
            raise ValueError("Implement _get_sparsedata for Matrix Datatype {0:s}".format(type(m)))

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
        self.var["S"] = self.nu * inner(grad(u), grad(w_test)) * dx
        self.var["K"] = inner(grad(self.u_stat) * u, w_test) * dx
        self.var["R"] = inner(grad(u) * self.u_stat, w_test) * dx
        self.var["G"] = p * div(w_test) * dx
        self.var["GT"] = div(u) * p_test * dx

        for (controlfunc, indices, name) in self.const.ASSEMBLER_BOUNDARY_CONTROLS:
            self.var["B_{0:s}".format(name)] = \
                self.nu*(1.0 / self.penalty_eps) * inner(controlfunc, w_test) * ds(indices)
                # (1.0 / self.penalty_eps) * inner(controlfunc, w_test) * ds(indices)

            self.var["M_{0:s}".format(name)] = \
                self.nu*(1.0 / self.penalty_eps) * inner(u, w_test) * ds(indices)
                # (1.0 / self.penalty_eps) * inner(u, w_test) * ds(indices)


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

        # for p in self.const.ASSEMBLER_OBSERVER_POINTS:
        #    print "P1({0:.2f}, {1:.2f})".format(p[0], p[1])

        points = [Point(*x) for x in self.const.ASSEMBLER_OBSERVER_POINTS]
        dists = dict([(p, float("inf")) for p in points])
        idxs = dict([(p, 0) for p in points])

        # find cell indices with midpoint has minimal distance
        for c in cells(self.mesh):
            for p in points:
                cdist = c.midpoint().distance(p)
                if cdist < dists[p]:
                    dists[p] = cdist
                    idxs[p] = c.index()

        # collect vertical dofs of cells with midpoint has minimal distance
        verticaldofs = np.empty(0, dtype=np.int64)
        for idx in idxs.values():
            verticaldofs = np.union1d(verticaldofs, self.V.sub(1).dofmap().cell_dofs(idx))

        # buid matrix
        C = scsp.lil_matrix((verticaldofs.size, self.V.dim()))
        j = 0
        for i in verticaldofs:
            C[j, i] = 1
            j += 1

        return C.tocsc()



    def _lns_npsc(self):
        if not self.ublas:
            raise ValueError("Call unparameterized_lns_ublas( to initialize attribute ublas.")

        self.npsc = {}
        for name, form in self.ublas.iteritems():
            if isinstance(form, dolfin.cpp.la.Vector):
                self.npsc[name] = form.array().reshape(form.array().size, 1)
            elif isinstance(form, dolfin.cpp.la.Matrix):
                self.npsc[name] = self._get_sparsedata(form)
            else:
                raise ValueError("Unknown type {0:s}".format(type(form)))

        # build input matrix B, concatenate or sum
        if self.const.NAME == "karman":
            self.npsc["B"] = np.concatenate([self.npsc[x] for x in self.npsc.keys() if x.startswith("B_")], axis=1)
        elif self.const.NAME == "bws":
            self.npsc["B"] = sum([self.npsc[x] for x in self.npsc.keys() if x.startswith("B_")])
        else:
            raise ValueError("unknown scneario name = {0:s}".format(self.const.NAME))

        # build matrix M of Boundary Controls
        self.npsc["M_BOUNDARY_CTRL"] = sum([self.npsc[x] for x in self.npsc.keys() if x.startswith("M_")])

        # build C
        self.npsc["C"] = self._buildC()

    def save_mtx(self):

        if not self.npsc:
            raise ValueError("Call unparameterized_lns_npsc to initialize attribute npsc.")

        for name, mat in self.npsc.iteritems():
            filename = self.const.ASSEMBLER_NAME_MTX(self.ref, name, self.RE)
            createdir(filename)
            #with open(filename, "w") as handle:
            #    if hasattr(mat, "eliminate_zeros"):
            #        mat.eliminate_zeros()
            #        scio.mmwrite(handle, mat, comment="{0:s},ref={1:d},RE={2:d}".format(name, self.ref, self.RE))
            if hasattr(mat, "eliminate_zeros"):
                mat.eliminate_zeros()
            write_matrix(filename, mat, "{0:s},ref={1:d},RE={2:d}".format(name, self.ref, self.RE))


    def save_mat(self):

        if not self.npsc:
            raise ValueError("Call unparameterized_lns_npsc to initialize attribute npsc.")
        filename = self.const.ASSEMBLER_MAT(self.ref, self.RE)
        createdir(filename)

        with open(filename, "w") as handle:
            scio.savemat(handle, self.npsc, do_compression=True)













