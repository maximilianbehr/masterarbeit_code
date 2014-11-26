from scipy.sparse.csr import csr_matrix
from scipy.sparse import lil_matrix

import numpy as np

from dolfin import *
from src.outputhandler import KarmanOutputHandler
from src.outputhandler import ProblemSolverOutputHandler
from src.problems.problem_mesh.karman import circle
from src.problems.problem_mesh.karman import GammaBallCtrlLower
from src.problems.problem_mesh.karman import GammaBallCtrlUpper


ref = 2
RE = 400
nu = 1.0 / RE
penalty_eps = 0.01

kohandler = KarmanOutputHandler()
psohandler = ProblemSolverOutputHandler("karman", "stat_newton")

OPTIONS = {
    "RE": 100,
    "mesh": kohandler.karman_mesh_xml(ref),
    "boundaryfunction": kohandler.karman_boundary_xml(ref),
    "u_stat": psohandler.u_xml(ref, RE),
    "penalty_eps": 0.01
}


# set matrix reordering off
parameters['reorder_dofs_serial'] = False


class Observer(SubDomain):
    """Domain for building submesh"""

    def inside(self, x, on_boundary):
        return between(x[0], (3.0, 3.5))



class LQR_Assembler():

    def __init__(self,argoptions):
        self.options = argoptions.copy
        self.mesh = Mesh(self.options["mesh"])
        self.boundaryfunction = MeshFunction("size_t", self.mesh, self.options["boundaryfunction"])
        self.g = Expression(("1/r * (x[0]-x0)", "1/r * (x[1]-y0)"), r=circle["r"], x0=circle["x0"], y0=circle["y0"])
        self.V = VectorFunctionSpace(self.mesh, "CG", 2)
        self.Q = FunctionSpace(self.mesh, "CG", 1)

        psohandler = ProblemSolverOutputHandler("karman", "stat_newton")
        self.u_stat = Function(self.V, psohandler.u_xml(ref, RE))


    def _get_sparsedata(self,m):
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
        dudt = TrialFunction(V)
        u = TrialFunction(V)
        p = TrialFunction(Q)

        # test functions
        w_test = TestFunction(V)
        p_test = TestFunction(Q)

        #weak formulation of linearized navier stokes eqn
        ds = Measure('ds')[boundaryfunction]
        self.var = {}

        self.var["M"] = inner(dudt, w_test) * dx
        self.var["S"] = inner(grad(u), grad(w_test)) * dx  #1/nu spaeter ranmultiplizieren
        self.var["B_upper"] = inner(g, w_test) * ds(GammaBallCtrlUpper.index)  #1/nu 1/eps spaeter ran
        self.var["B_lower"] = inner(g, w_test) * ds(GammaBallCtrlLower.index)  #1/nu 1/eps spaeter ran
        self.var["M_upper"] = inner(u, w_test) * ds(GammaBallCtrlUpper.index)  #1/nu 1/eps spaeter ran
        self.var["M_lower"] = inner(u, w_test) * ds(GammaBallCtrlLower.index)  #1/nu 1/eps spaeter ran
        self.var["K"] = inner(grad(u_stat) * u, w_test) * dx
        self.var["R"] = inner(grad(u) * u_stat, w_test) * dx
        self.var["G"] = -1 * p * div(w_test) * dx

        #nebenbedingung mit -1 multiplizieren
        self.var["Gt"] = -1 * div(u) * p_test * dx

    def unparameterized_lns_petsc(self):

        if not self.var:
            raise ValueError("Call unparameterized_lns_variational to initialize attribute var.")

        #store backend
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

        #store backend
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
            if key == "B_upper" or key == "B_lower":
                self.npsc[key] = val.array().reshape(val.array().size, 1)
            else:
                self.npsc[key] = self._get_sparsedata(val)

        self.npsc["B"] = np.concatenate((self.npsc["B_upper"], self.npsc["B_lower"]), axis=1)

        #build matrix C
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

        self.npsc["C"]=C.todense()



from pycmess import options
from pycmess import equation_dae2
from pycmess import PYCMESS_OP_TRANSPOSE
from pycmess import lrnm
from pycmess import lrnm_res

# create opt instance
opt = options()
opt.adi.output = 1;
opt.nm.output = 1;
opt.nm.res2_tol = 1e-5;

eqn = equation_dae2()
eqn.M = mat_numpyscipy["M"]
eqn.A = 1.0 / nu * (
-mat_numpyscipy["S"] - 1.0 / penalty_eps * (mat_numpyscipy["M_lower"] + mat_numpyscipy["M_upper"])) - mat_numpyscipy[
    "K"] - mat_numpyscipy["R"]
eqn.G = -1 * mat_numpyscipy["G"]
eqn.B = 1.0 / nu * 1.0 / penalty_eps * mat_numpyscipy["B"]
eqn.C = Ccsr
eqn.delta = -0.02

#solve generalized riccati equation
opt.adi.type = PYCMESS_OP_TRANSPOSE
result = lrnm(eqn, opt)
Z = result[0]
res2 = result[1]
iter = result[2]
print res2[res2.size - 1]
result = lrnm_res(eqn, opt, Z)
res = result[0]
rel = result[1]
print "rel= %e \t res= %e \t res/rel= %e \n" % (rel, res, res / rel)