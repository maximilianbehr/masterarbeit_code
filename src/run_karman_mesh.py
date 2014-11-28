from dolfin import parameters
from dolfin.cpp.common import begin
from dolfin.cpp.common import end
from dolfin.cpp.common import info
from problems.problem_mesh.karman import MeshBuilder
from outputhandler import KarmanOutputHandler


refalgs = ["bisection", "iterative_bisection", "recursive_bisection", "regular_cut"]
refs = range(1, 8)

if __name__ == "__main__":

    # object for saving mesh and boundary function
    kohandler = KarmanOutputHandler()

    for refalg in refalgs:

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

            #save mesh
            kohandler.karman_save_mesh(meshbuilder)

            #save boundaryfunction
            kohandler.karman_save_boundaryfunction(meshbuilder)

        end()



