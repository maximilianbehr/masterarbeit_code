from dolfin import *



# parameters["reorder_dofs_serial"]=False



mesh = UnitSquareMesh(10, 10)


class Observer(SubDomain):
    """Left Boundary"""

    def inside(self, x, on_boundary):
        return between(x[0], (1. / 4, 3. / 4))


ob = Observer()
submesh = SubMesh(mesh, ob)

indices_cells_in_observer = submesh.data().array("parent_cell_indices", 2)

V = VectorFunctionSpace(mesh, "CG", 2)

for i in indices_cells_in_observer:
    print V.sub(0).dofmap().cell_dofs(i), V.sub(1).dofmap().cell_dofs(i), V.dofmap().cell_dofs(i)

