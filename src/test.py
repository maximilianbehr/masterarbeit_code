from dolfin import *
from scipy.sparse import csr_matrix
from scipy.sparse import lil_matrix
import numpy as np

def tocsr(m):
    rowptr, colptr, vals = m.data()
    rows = m.size(0)
    cols = m.size(1)
    return csr_matrix((vals, colptr, rowptr), (rows, cols))


parameters["linear_algebra_backend"] = "uBLAS"

mesh = Mesh('/Users/daniels/PycharmProjects/master/data/1.3.0+/karman/recursive_bisection/ref_2/mesh.xml.gz')

V = VectorFunctionSpace(mesh, "CG", 2)

p1 = Point(3.5, 0.25)
p2 = Point(3.5, 0.75)
points = [p1, p2]

dists = {p1: float("inf"), p2: float("inf")}
idxs = {p1: 0, p2: 0}

# find cell indices with midpoint has minimal distance
for c in cells(mesh):
    for p in points:
        cdist = c.midpoint().distance(p)
        if cdist < dists[p]:
            dists[p] = cdist
            idxs[p] = c.index()

#collect vertical dofs of cells with midpoint has minimal distance
verticaldofs = np.empty(0, dtype=np.int32)
for idx in idxs.values():
    verticaldofs = np.union1d(verticaldofs, V.sub(1).dofmap().cell_dofs(idx))

C = lil_matrix((verticaldofs.size, V.dim()))
j = 0
for i in verticaldofs:
    C[j, i] = 1
    j += 1

print C



#if indices.size == 0:
#    raise ValueError("There are no cells with dofs in the domain. Change inside Method in Observer")

#C = lil_matrix((indices.size, self.V.dim()))
#j = 0
#for i in indices:
#    C[j, i] = 1
#    j += 1



# dof reordering off
# parameters["reorder_dofs_serial"] = False
# V = VectorFunctionSpace(mesh, "CG", 2)
#Q = FunctionSpace(mesh, "CG", 1)
#q = TestFunction(Q)
#p = TrialFunction(Q)
#Qd = Q.dofmap()
#M_var = inner(p, q) * dx
#M = assemble(M_var)
#M_csr = tocsr(M)


# dof reordering on

# parameters["reorder_dofs_serial"] = True
#V2 = VectorFunctionSpace(mesh, "CG", 2)
#Q2 = FunctionSpace(mesh, "CG", 1)
#q2 = TestFunction(Q2)
#p2 = TrialFunction(Q2)
#Qd2 = Q2.dofmap()
#M2_var = inner(p2, q2) * dx
#M2 = assemble(M2_var)
#M2_csr = tocsr(M2)

