from dolfin import *
from scipy.sparse import csr_matrix


def tocsr(m):
    rowptr, colptr, vals = m.data()
    rows = m.size(0)
    cols = m.size(1)
    return csr_matrix((vals, colptr, rowptr), (rows, cols))


parameters["linear_algebra_backend"] = "uBLAS"

mesh = Mesh('/Users/daniels/PycharmProjects/master/data/1.3.0+/karman/recursive_bisection/macro/mesh.xml.gz')

# dof reordering off
parameters["reorder_dofs_serial"] = False
V = VectorFunctionSpace(mesh, "CG", 2)
Q = FunctionSpace(mesh, "CG", 1)
q = TestFunction(Q)
p = TrialFunction(Q)
Qd = Q.dofmap()
M_var = inner(p, q) * dx
M = assemble(M_var)
M_csr = tocsr(M)


# dof reordering on

parameters["reorder_dofs_serial"] = True
V2 = VectorFunctionSpace(mesh, "CG", 2)
Q2 = FunctionSpace(mesh, "CG", 1)
q2 = TestFunction(Q2)
p2 = TrialFunction(Q2)
Qd2 = Q2.dofmap()
M2_var = inner(p2, q2) * dx
M2 = assemble(M2_var)
M2_csr = tocsr(M2)

