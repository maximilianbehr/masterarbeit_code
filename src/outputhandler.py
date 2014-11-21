import os
import socket
from dolfin import __version__
from dolfin.cpp.io import File
from dolfin import parameters


extharddrive = "/media/UNTITLED/"
extharddrivemac = "/Volumes/UNTITLED/"


class KarmanOutputHandler():

    def karman_outputdir(self):
        """return outputdir depending on pc name and existence of external harddrive
        for karman meshes and boundary functions
        """
        refinementalg = parameters["refinement_algorithm"]

        dirname = os.path.join("data",__version__,"karman",refinementalg)

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
            if os.path.isdir(extharddrivemac):
                return os.path.join(extharddrivemac, dirname)
            else:
                return os.path.abspath(os.path.join("/Users/daniels/PycharmProjects/master/", dirname))
        else:
            raise NotImplementedError("Output Path for %s not implemented" % hostname)


    def _karman_file(self, name, num):
        outputdir = self.karman_outputdir()
        if num:
            return os.path.join(outputdir, "ref_%s" % num, name)
        return os.path.join(outputdir, "macro", name)

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
        File(self.karman_mesh_xml(meshbuilder.refinelevel),"compressed") << meshbuilder.mesh
        File(self.karman_mesh_pvd(meshbuilder.refinelevel),"compressed") << meshbuilder.mesh
        File(self.karman_mesh_xdmf(meshbuilder.refinelevel) ,"compressed")<< meshbuilder.mesh

    def karman_save_boundaryfunction(self, meshbuilder):
        File(self.karman_boundary_xdmf(meshbuilder.refinelevel),"compressed") << meshbuilder.boundaryfunction
        File(self.karman_boundary_xml(meshbuilder.refinelevel),"compressed") << meshbuilder.boundaryfunction
        File(self.karman_boundary_pvd(meshbuilder.refinelevel),"compressed") << meshbuilder.boundaryfunction


class ProblemSolverOutputHandler():
    def __init__(self, problemname, solvername):
        self.problemname = problemname
        self.solvername = solvername

    def outputdir(self):
        """return outputlocation depending on pc name and existence of external
        for a specific problem and solver"""


        refinementalg = parameters["refinement_algorithm"]
        dirname = os.path.join("results",__version__,self.problemname,self.solvername,refinementalg)

        hostname = socket.gethostname()
        if hostname == "pc747":
            if os.path.isdir(extharddrive):
                return os.path.join(extharddrive, dirname)
            return os.path.abspath(os.path.join("/scratch/behr/master/src/", dirname))

        elif hostname == "pc800":
            if os.path.isdir(extharddrive):
                return os.path.join(extharddrive, dirname)
            return os.path.abspath(os.path.join("/scratch/behr/master/src/", dirname))

        elif hostname == "pc633":
            if os.path.isdir(extharddrive):
                return os.path.join(extharddrive, dirname)
            return os.path.abspath(os.path.join("/scratch/behr/master/src/", dirname))

        elif hostname == "jack":
            if os.path.isdir(extharddrivemac):
                return os.path.join(extharddrivemac, dirname)
            return os.path.abspath(os.path.join("/Users/daniels/PycharmProjects/master/", dirname))

        else:
            raise NotImplementedError("Output Path for %s not implemented" % hostname)


    def _file(self, name, num, RE):
        assert(num > 0)
        outputdir = self.outputdir()
        return os.path.join(outputdir, "ref_%s" % num, "RE_{0}".format(RE), name)

    def u_xml(self, num, RE):
        name = "u.xml"
        return self._file(name, num, RE)

    def u_t_xml(self, num, RE):
        name = "u_{0:1.2e}.xml"
        return self._file(name, num,RE)

    def p_xml(self, num, RE):
        name = "p.xml"
        return self._file(name, num, RE)

    def p_t_xml(self, num, RE):
        name = "u_{0:1.2e}.xml"
        return self._file(name, num,RE)

    def u_pvd(self, num, RE):
        name = "u.pvd"
        return self._file(name, num, RE)

    def p_pvd(self, num, RE):
        name = "p.pvd"
        return self._file(name, num, RE)

    def options_json(self,num,RE):
        name = "options.json"
        return self._file(name, num, RE)