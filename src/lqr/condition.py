import scipy.io as scio
import scipy.sparse as scsp
import scipy.sparse.linalg as scspli
import scipy.linalg as scla
import numpy as np
from pyamg.util.linalg import condest

# set floating point error handling
#np.seterr(all="raise", divide="raise", over="raise", under="raise", invalid="raise")
np.seterr(all="ignore")



class Condition():
    def __init__(self, const, ref, RE):

        # parameters
        self.ref = ref
        self.RE = RE
        self.const = const

        # read compress system
        self.names = ["M", "M_BOUNDARY_CTRL", "S", "R", "K", "G", "GT"]
        self.ctrlmat = {}
        self.simmat = {}
        for name in self.names:
            self.simmat[name] = scio.mmread(const.ASSEMBLER_COMPRESS_SIM_NAME_MTX(ref, name, RE))
            self.ctrlmat[name] = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_NAME_MTX(ref, name, RE))

        # build system matrices sim
        self.nv_sim = self.simmat["G"].shape[0]
        self.np_sim = self.simmat["G"].shape[1]
        self.Asys_sim = -self.simmat["S"]-self.simmat["R"]-self.simmat["K"]-self.simmat["M_BOUNDARY_CTRL"]
        u = scsp.hstack([self.Asys_sim, self.simmat["G"]])
        l = scsp.hstack([self.simmat["GT"], scsp.csr_matrix((self.np_sim, self.np_sim))])
        self.Asim = scsp.vstack([u, l])

        # build system matrices ctrl
        self.nv_ctrl = self.ctrlmat["G"].shape[0]
        self.np_ctrl = self.ctrlmat["G"].shape[1]
        self.Asys_ctrl = -self.ctrlmat["S"]-self.ctrlmat["R"]-self.ctrlmat["K"]-self.ctrlmat["M_BOUNDARY_CTRL"]
        u = scsp.hstack([self.Asys_ctrl, self.ctrlmat["G"]])
        l = scsp.hstack([self.ctrlmat["GT"], scsp.csr_matrix((self.np_ctrl, self.np_ctrl))])
        self.Actrl = scsp.vstack([u, l])

    def compute_condests(self):
        print "ref={0:d} RE={1:d}".format(self.ref, self.RE)

        for name in self.names:
            if name != "G" and name != "GT":
                print "Simulate cond({0:s})\t={1:e}".format(name, condest(self.simmat[name]))
        print "Simulate cond(Asys_sim)\t={1:e}".format(name, condest(self.Asys_sim))
        print "Simulate cond(Asim)\t={1:e}".format(name, condest(self.Asim))

        print "------------------------------------------------------------------------"
        for name in self.names:
            if name != "G" and name != "GT":
                print "Control cond({0:s})\t={1:e}".format(name, condest(self.ctrlmat[name]))
        print "Simulate cond(Asys_ctrl)\t={1:e}".format(name, condest(self.Asys_ctrl))
        print "Control cond(Actrl)\t={1:e}".format(name, condest(self.Actrl))



