import os

from scipy.sparse.csr import csr_matrix
from scipy.sparse import lil_matrix
from scipy.io import mmwrite
from scipy.io import savemat
import numpy as np
from dolfin import *
from src.problems.problem_mesh.karman import circle
from src.problems.problem_mesh.karman import GammaBallCtrlLower
from src.problems.problem_mesh.karman import GammaBallCtrlUpper




# set matrix reordering off
parameters['reorder_dofs_serial'] = False


class Observer(SubDomain):
    """Domain for building submesh"""

    def inside(self, x, on_boundary):
        return between(x[0], (3.0, 3.5))


class LQR_Assembler():
    def __init__(self, argoptions):
        self.options = argoptions.copy()
        self.mesh = Mesh(self.options["mesh"])
        self.boundaryfunction = MeshFunction("size_t", self.mesh, self.options["boundaryfunction"])
        self.g = Expression(("1/r * (x[0]-x0)", "1/r * (x[1]-y0)"), r=circle["r"], x0=circle["x0"], y0=circle["y0"])
        self.V = VectorFunctionSpace(self.mesh, "CG", 2)
        self.Q = FunctionSpace(self.mesh, "CG", 1)

        self.u_stat = Function(self.V, self.options["u_stat"])


    def _get_sparsedata(self, m):
        """Function returns a scipy CSR Matrix for given Matrix (dolfin) M"""
        try:
            rowptr, colptr, vals = m.data()
            rows = m.size(0)
            cols = m.size(1)
            return csr_matrix((vals, colptr, rowptr), (rows, cols))
        except:
            raise ValueError(u'Implement _get_sparsedata for Matrix Datatype {0:s}'.format(type(m)))


    def unparameterized_lns_variational(self):
        # trial functions
        dudt = TrialFunction(self.V)
        u = TrialFunction(self.V)
        p = TrialFunction(self.Q)

        # test functions
        w_test = TestFunction(self.V)
        p_test = TestFunction(self.Q)

        # weak formulation of linearized navier stokes eqn
        ds = Measure('ds')[self.boundaryfunction]
        self.var = {}

        self.var["M"] = inner(dudt, w_test) * dx
        self.var["S"] = self.options["RE"] * inner(grad(u), grad(w_test)) * dx
        self.var["Bupper"] = self.options["RE"] * inner(self.g, w_test) * ds(
            GammaBallCtrlUpper.index)  #1/eps spaeter ran
        self.var["Blower"] = self.options["RE"] * inner(self.g, w_test) * ds(
            GammaBallCtrlLower.index)  #1/eps spaeter ran
        self.var["Mupper"] = self.options["RE"] * inner(u, w_test) * ds(GammaBallCtrlUpper.index)  #1/eps spaeter ran
        self.var["Mlower"] = self.options["RE"] * inner(u, w_test) * ds(GammaBallCtrlLower.index)  #1/eps spaeter ran
        self.var["K"] = inner(grad(self.u_stat) * u, w_test) * dx
        self.var["R"] = inner(grad(u) * self.u_stat, w_test) * dx
        self.var["G"] = -1 * p * div(w_test) * dx

        #nebenbedingung mit -1 multiplizieren
        self.var["Gt"] = -1 * div(u) * p_test * dx

    def unparameterized_lns_petsc(self):

        if not self.var:
            raise ValueError("Call unparameterized_lns_variational to initialize attribute var.")

        # store backend
        la_backend = parameters["linear_algebra_backend"]

        # assemble matrices for PETSC backend
        parameters["linear_algebra_backend"] = "PETSc"
        self.petsc = {}
        for key, val in self.var.iteritems():
            self.petsc[key] = assemble(val)

        #restore backend
        parameters["linear_algebra_backend"] = la_backend

    def unparameterized_lns_ublas(self):
        if not self.var:
            raise ValueError("Call unparameterized_lns_variational to initialize attribute var.")

        # store backend
        la_backend = parameters["linear_algebra_backend"]

        # assemble matrices for PETSC backend
        parameters["linear_algebra_backend"] = "uBLAS"
        self.ublas = {}
        for key, val in self.var.iteritems():
            self.ublas[key] = assemble(val)

        #restore backend
        parameters["linear_algebra_backend"] = la_backend

    def unparameterized_lns_npsc(self):
        if not self.ublas:
            raise ValueError("Call unparameterized_lns_ublas( to initialize attribute ublas.")

        self.npsc = {}
        for key, val in self.ublas.iteritems():
            if key == "Bupper" or key == "Blower":
                self.npsc[key] = val.array().reshape(val.array().size, 1)
            else:
                self.npsc[key] = self._get_sparsedata(val)

        self.npsc["B"] = np.concatenate((self.npsc["Bupper"], self.npsc["Blower"]), axis=1)

        # build matrix C
        ob = Observer()
        submesh = SubMesh(self.mesh, ob)
        indices = np.empty(0, dtype=np.int32)
        for i in submesh.data().array("parent_cell_indices", 2):
            indices = np.union1d(indices, self.V.sub(1).dofmap().cell_dofs(i))

        if indices.size == 0:
            raise ValueError("There are no cells with dofs in the domain. Change inside Method in Observer")

        C = lil_matrix((indices.size, self.V.dim()))
        j = 0
        for i in indices:
            C[j, i] = 1
            j += 1

        self.npsc["C"] = C.todense()


    def unparameterized_lns_mtx(self):

        if not self.npsc:
            raise ValueError("Call unparameterized_lns_npsc to initialize attribute npsc.")

        if not os.path.exists(os.path.dirname(self.options["M_mtx"])):
            os.makedirs(os.path.dirname(self.options["M_mtx"]))

        with open(self.options["M_mtx"], "w") as handle:
            mmwrite(handle, self.npsc["M"])

        with open(self.options["S_mtx"], "w") as handle:
            mmwrite(handle, self.npsc["S"])

        with open(self.options["Mlower_mtx"], "w") as handle:
            mmwrite(handle, self.npsc["Mlower"])

        with open(self.options["Mupper_mtx"], "w") as handle:
            mmwrite(handle, self.npsc["Mupper"])

        with open(self.options["K_mtx"], "w") as handle:
            mmwrite(handle, self.npsc["K"])

        with open(self.options["R_mtx"], "w") as handle:
            mmwrite(handle, self.npsc["R"])

        with open(self.options["G_mtx"], "w") as handle:
            mmwrite(handle, self.npsc["G"])

        with open(self.options["Gt_mtx"], "w") as handle:
            mmwrite(handle, self.npsc["Gt"])

        with open(self.options["Blower_mtx"], "w") as handle:
            mmwrite(handle, self.npsc["Blower"])

        with open(self.options["Bupper_mtx"], "w") as handle:
            mmwrite(handle, self.npsc["Bupper"])

        with open(self.options["B_mtx"], "w") as handle:
            mmwrite(handle, self.npsc["B"])

        with open(self.options["C_mtx"], "w") as handle:
            mmwrite(handle, self.npsc["C"])


    def unparametrized_lns_mat(self):

        if not self.npsc:
            raise ValueError("Call unparameterized_lns_npsc to initialize attribute npsc.")

        with open(self.options["mat"], "w") as handle:
            savemat(handle, self.npsc, do_compression=True)

    def save_options(self):
        try:
            import json

            fname = self.options["options_json"]
            if not os.path.exists(os.path.dirname(fname)):
                os.makedirs(os.path.dirname(fname))
            with open(fname, "w") as handle:
                json.dump(self.options, handle)
        except ImportError:
            print "Cannot import json. options not stored"












