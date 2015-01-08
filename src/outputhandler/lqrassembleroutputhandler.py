
from src.outputhandler import extharddrive
from src.outputhandler import extharddrivemac

import os
import socket
import sys

from dolfin import __version__
from dolfin import parameters

class LQRAssemblerOutputHandler():
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

        dirname = os.path.join("results", __version__, "karman", "lqr_assemble", \
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

    def M_mtx(self):
        name = "M.mtx"
        return self._file(name)

    def S_mtx(self):
        name = "S.mtx"
        return self._file(name)

    def Mlower_mtx(self):
        name = "M_lower.mtx"
        return self._file(name)

    def Mupper_mtx(self):
        name = "M_upper.mtx"
        return self._file(name)

    def K_mtx(self):
        name = "K.mtx"
        return self._file(name)

    def R_mtx(self):
        name = "R.mtx"
        return self._file(name)

    def G_mtx(self):
        name = "G.mtx"
        return self._file(name)

    def Gt_mtx(self):
        name = "Gt.mtx"
        return self._file(name)

    def Blower_mtx(self):
        name = "B_lower.mtx"
        return self._file(name)

    def Bupper_mtx(self):
        name = "B_upper.mtx"
        return self._file(name)

    def B_mtx(self):
        name = "B.mtx"
        return self._file(name)

    def C_mtx(self):
        name = "C.mtx"
        return self._file(name)

    def mat(self):
        name = "lns.mat"
        return self._file(name)

    def options_json_assembler(self):
        name = "options_assembler.json"
        return self._file(name)

    def options_json_solver(self):
        name = "options_solver.json"
        return self._file(name)

    def log_assembler(self):
        name = "log_assembler.log"
        return self._file(name)

    def log_solver(self):
        name = "log_solver.log"
        return self._file(name)

    def Z_mtx(self):
        name = "Z.mtx"
        return self._file(name)

    def Kinf_mtx(self):
        name = "Kinf.mtx"
        return self._file(name)

    def Feed0_mtx(self):
        name = "Feed0.mtx"
        return self._file(name)

    def Feed1_mtx(self):
        name = "Feed1.mtx"
        return self._file(name)

    def res2_txt(self):
        name = "lqr_res2.txt"
        return self._file(name)

    def M_PETSC(self):
        name = "M.PETSC"
        return self._file(name)

    def S_PETSC(self):
        name = "S.PETSC"
        return self._file(name)

    def Mlower_PETSC(self):
        name = "M_lower.PETSC"
        return self._file(name)

    def Mupper_PETSC(self):
        name = "M_upper.PETSC"
        return self._file(name)

    def K_PETSC(self):
        name = "K.PETSC"
        return self._file(name)

    def R_PETSC(self):
        name = "R.PETSC"
        return self._file(name)

    def G_PETSC(self):
        name = "G.PETSC"
        return self._file(name)

    def Gt_PETSC(self):
        name = "Gt.PETSC"
        return self._file(name)

    def Blower_PETSC(self):
        name = "B_lower.PETSC"
        return self._file(name)

    def Bupper_PETSC(self):
        name = "B_upper.PETSC"
        return self._file(name)

    def B_PETSC(self):
        name = "B.PETSC"
        return self._file(name)

    def C_PETSC(self):
        name = "C.PETSC"
        return self._file(name)