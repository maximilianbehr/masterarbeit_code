from src.outputhandler import extharddrive
from src.outputhandler import extharddrivemac

import os
import socket
import sys

from dolfin import __version__
from dolfin import parameters


class ProblemSolverOutputHandler():
    def __init__(self, problemname, solvername, num, RE):
        assert (num >= 0)
        self.num = num
        self.RE = RE
        self.problemname = problemname
        self.solvername = solvername
        self.outputdir = self._outputdir()

    def _outputdir(self):
        """return outputlocation depending on pc name and existence of external
        for a specific problem and solver"""

        refinementalg = parameters["refinement_algorithm"]
        dirname = os.path.join("results", __version__, self.problemname, self.solvername, refinementalg, \
                               "ref_{0:d}".format(self.num), "RE_{0:d}".format(self.RE))

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
            if sys.platform=="linux2" or sys.platform=="linux":
                return os.path.abspath(os.path.join("/home/daniels/PycharmProjects/master/", dirname))
            else:
                if os.path.isdir(extharddrivemac):
                    return os.path.join(extharddrivemac, dirname)
                else:
                    return os.path.abspath(os.path.join("/Users/daniels/PycharmProjects/master/", dirname))

        else:
            raise NotImplementedError("Output Path for %s not implemented" % hostname)


    def _file(self, name, ):
        return os.path.join(self.outputdir, name)

    def u_xml(self):
        name = "u.xml"
        return self._file(name)

    def u_t_xml(self):
        name = "u_{0:1.2e}.xml"
        return self._file(name)

    def p_xml(self):
        name = "p.xml"
        return self._file(name)

    def p_t_xml(self):
        name = "p_{0:1.2e}.xml"
        return self._file(name)

    def u_pvd(self):
        name = "u.pvd"
        return self._file(name)

    def p_pvd(self):
        name = "p.pvd"
        return self._file(name)

    def options_json(self):
        name = "options.json"
        return self._file(name)

    def log(self):
        name = "log.log"
        return self._file(name)
