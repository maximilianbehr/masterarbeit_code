import os
import socket

from dolfin import parameters
from dolfin import __version__
from dolfin.cpp.io import File
from dolfin.cpp.common import begin
from dolfin.cpp.common import end
from dolfin.cpp.common import info
from problems.problem_mesh.karman import MeshBuilder

import ipdb

class MeshPathHandler():
    def __init__(self):
        self.extraharddrive = "/media/UNTITLED/"
        self.extraharddrivemac = "/Volumes/UNTITLED/"

    def getOutputlocation(self):
        """return outputlocation depending on pc name and existence of external harddrive
        """

        hostname = socket.gethostname()
        if hostname == "pc747":
            if os.path.isdir(self.extraharddrive):
                return os.path.join(self.extraharddrive, "data/karman")
            else:
                return os.path.abspath("/scratch/behr/master/src/data/karman")
        elif hostname == "pc800":
            if os.path.isdir(self.extraharddrive):
                return os.path.join(self.extraharddrive, "data/karman")
            else:
                return os.path.abspath("/scratch/behr/master/src/data/karman")
        elif hostname == "pc633":
            if os.path.isdir(self.extraharddrive):
                return os.path.join(self.extraharddrive, "data/karman")
            else:
                return os.path.abspath("/scratch/behr/master/src/data/karman")
        elif hostname == "jack":
            if os.path.isdir(self.extraharddrivemac):
                return os.path.join(self.extraharddrivemac, "data/karman")
            else:
                return os.path.abspath("/Users/daniels/PycharmProjects/master/data/karman")

    def _karman_file(self, name, num):
        outputlocation = self.getOutputlocation()
        refinementalg = parameters["refinement_algorithm"]
        if num:
            return File(os.path.join(outputlocation, __version__, refinementalg, "ref_%s" % num, name), "compressed")
        return File(os.path.join(outputlocation, __version__, "macro", name), "compressed")

    def karman_mesh_xml(self, num):
        name = "mesh.xml.gz"
        return self._karman_file(name, num)

    def karman_boundary_xml(self, num):
        name = "boundary.xml.gz"
        return self._karman_file(name, num)

    def karman_mesh_xdmf(self, num):
        name = "problem_mesh.xdmf"
        return self._karman_file(name, num)

    def karman_boundary_xdmf(self, num):
        name = "boundary.xdmf"
        return self._karman_file(name, num)

    def karman_mesh_pvd(self, num):
        name = "problem_mesh.pvd"
        return self._karman_file(name, num)

    def karman_boundary_pvd(self, num):
        name = "boundary.pvd"
        return self._karman_file(name, num)

    def save_mesh(self, meshbuilder):
        self.karman_mesh_xml(meshbuilder.refinelevel) << meshbuilder.mesh
        self.karman_mesh_pvd(meshbuilder.refinelevel) << meshbuilder.mesh
        self.karman_mesh_xdmf(meshbuilder.refinelevel) << meshbuilder.mesh

    def save_boundaryfunction(self, meshbuilder):
        self.karman_boundary_xdmf(meshbuilder.refinelevel) << meshbuilder.boundaryfunction
        self.karman_boundary_xml(meshbuilder.refinelevel) << meshbuilder.boundaryfunction
        self.karman_mesh_pvd(meshbuilder.refinelevel) << meshbuilder.boundaryfunction


if __name__ == "__main__":

    # object for saving mesh and boundary function
    meshpathhandler = MeshPathHandler()

    for refalg in ["bisection", "iterative_bisection", "recursive_bisection", "regular_cut"]:

        parameters["refinement_algorithm"] = refalg

        # build macro mesh
        meshbuilder = MeshBuilder()

        #save mesh and meshfunction
        meshpathhandler.save_mesh(meshbuilder)
        meshpathhandler.save_boundaryfunction(meshbuilder)

        begin("Refinement with %s" % parameters["refinement_algorithm"])
        for ref in range(1, 5):
            info("Level %d" % ref)

            #refine mesh
            meshbuilder.refine()

            #save mesh
            meshpathhandler.save_mesh(meshbuilder)

            #save boundaryfunction
            meshpathhandler.save_boundaryfunction(meshbuilder)

        end()



