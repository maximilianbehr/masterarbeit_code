from dolfin import parameters
from dolfin.cpp.common import begin
from dolfin.cpp.common import end
from dolfin.cpp.common import info
from problems.problem_mesh.karman import MeshBuilder
from src.outputhandler.karmanoutputhandler import KarmanOutputHandler

# set dof reordering off
parameters["reorder_dofs_serial"] = False
parameters["form_compiler"]["optimize"] = True

# refalgs = ["bisection", "iterative_bisection", "recursive_bisection", "regular_cut"]
refalgs = ["recursive_bisection"]

# refs = range(1, 7)
refs = range(1, 5)

if __name__ == "__main__":

    # object for saving mesh and boundary function
    kohandler = KarmanOutputHandler()

    for refalg in refalgs:

        # set refinement algorithm
        parameters["refinement_algorithm"] = refalg

        # build macro mesh
        meshbuilder = MeshBuilder()

        # save macro mesh and meshfunction
        kohandler.karman_save_mesh(meshbuilder)
        kohandler.karman_save_boundaryfunction(meshbuilder)

        begin("Refinement with %s" % parameters["refinement_algorithm"])
        for ref in refs:
            info("Refinement {0}".format(ref))

            # refine mesh
            meshbuilder.refine()

            # save mesh
            kohandler.karman_save_mesh(meshbuilder)

            # save boundaryfunction
            kohandler.karman_save_boundaryfunction(meshbuilder)

        end()


