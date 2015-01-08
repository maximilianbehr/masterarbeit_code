from src.outputhandler import extharddrive
from src.outputhandler import extharddrivemac

import os
import socket
import sys

from dolfin import __version__
from dolfin.cpp.io import File
from dolfin import parameters


class KarmanOutputHandler():
    def __init__(self):
        self.outputdir = self._karman_outputdir()

    def _karman_outputdir(self):
        """return outputdir depending on pc name and existence of external harddrive
        for karman meshes and boundary functions
        """
        refinementalg = parameters["refinement_algorithm"]
        dirname = os.path.join("data", __version__, "karman", refinementalg)

        hostname = socket.gethostname()
        if hostname == "pc747":
            if os.path.isdir(extharddrive):
                return os.path.join(extharddrive, dirname)
            else:
                return os.path.abspath(os.path.join("/scratch/behr/master/src/", dirname))
        elif hostname == "pc800":
            if os.path.isdir(extharddrive):
                return os.path.join(extharddrive, dirname)
            else:
                return os.path.abspath(os.path.join("/scratch/behr/master/src/", dirname))
        elif hostname == "pc633":
            if os.path.isdir(extharddrive):
                return os.path.join(extharddrive, dirname)
            else:
                return os.path.abspath(os.path.join("/scratch/behr/master/src/", dirname))
        elif hostname == "jack":
            if sys.platform=="linux2" or sys.platform=="linux":
                return os.path.abspath(os.path.join("/home/daniels/PycharmProjects/master/", dirname))
            else:
                if os.path.isdir(extharddrivemac):
                    return os.path.join(extharddrivemac, dirname)
                else:
                    return os.path.abspath(os.path.join("/Users/daniels/PycharmProjects/master/", dirname))
        else:
            raise NotImplementedError("Output Path for %s not implemented" % hostname)

    def _karman_file(self, name, num):
        if num:
            return os.path.join(self.outputdir, "ref_%s" % num, name)
        return os.path.join(self.outputdir, "macro", name)

    def karman_mesh_xml(self, num):
        name = "mesh.xml.gz"
        return self._karman_file(name, num)

    def karman_boundary_xml(self, num):
        name = "boundary.xml.gz"
        return self._karman_file(name, num)

    def karman_mesh_xdmf(self, num):
        name = "mesh.xdmf"
        return self._karman_file(name, num)

    def karman_boundary_xdmf(self, num):
        name = "boundary.xdmf"
        return self._karman_file(name, num)

    def karman_mesh_pvd(self, num):
        name = "mesh.pvd"
        return self._karman_file(name, num)

    def karman_boundary_pvd(self, num):
        name = "boundary.pvd"
        return self._karman_file(name, num)

    def karman_save_mesh(self, meshbuilder):
        File(self.karman_mesh_xml(meshbuilder.refinelevel), "compressed") << meshbuilder.mesh
        File(self.karman_mesh_pvd(meshbuilder.refinelevel), "compressed") << meshbuilder.mesh
        File(self.karman_mesh_xdmf(meshbuilder.refinelevel), "compressed") << meshbuilder.mesh

    def karman_save_boundaryfunction(self, meshbuilder):
        File(self.karman_boundary_xdmf(meshbuilder.refinelevel), "compressed") << meshbuilder.boundaryfunction
        File(self.karman_boundary_xml(meshbuilder.refinelevel), "compressed") << meshbuilder.boundaryfunction
        File(self.karman_boundary_pvd(meshbuilder.refinelevel), "compressed") << meshbuilder.boundaryfunction
