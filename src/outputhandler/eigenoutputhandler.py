from src.outputhandler import extharddrive
from src.outputhandler import extharddrivemac

import os
import socket
import sys

from dolfin import __version__
from dolfin import parameters



class EigenOutputHandler():
    def __init__(self, num, RE):
        assert (num > 0)
        self.num = num
        self.RE = RE
        self.outputdir = self._outputdir()

    def _outputdir(self):
        """return outputdir depending on pc name and existence of external harddrive
        for karman meshes and boundary functions
        """
        refinementalg = parameters["refinement_algorithm"]

        dirname = os.path.join("results", __version__, "karman", "eigen", \
                               refinementalg, "ref_{0:d}".format(self.num), "RE_{0:d}".format(self.RE))

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

    def _file(self, name):
        return os.path.join(self.outputdir, name)

    def eig_eps(self):
        name = "eig.eps"
        return self._file(name)

    def eig_nopenalty_eps(self):
        name = "eig_nopenalty.eps"
        return self._file(name)

    def eig_bernoulli_eps(self):
        name = "eig_bernoulli.eps"
        return self._file(name)

    def eig_mtx(self):
        name = "eig.mtx"
        return self._file(name)

    def eig_nopenalty_mtx(self):
        name = "eig_nopenalty.mtx"
        return self._file(name)

    def eig_bernoulli_mtx(self):
        name = "eig_bernoulli.mtx"
        return self._file(name)

    def options_json(self):
        name = "options.json"
        return self._file(name)

    def log(self):
        name = "log.log"
        return self._file(name)