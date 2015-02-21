from linearized_aux import compress_mat
from linearized_aux import inner_outer_nodes
import scipy.io as scio
from dolfin import *
from src.aux import createdir, write_matrix
import numpy as np

class CompressAssembler():
    def __init__(self, const, ref, RE, flag):

        # parameters
        self.ref = ref
        self.RE = RE
        self.flag = flag
        self.const = const

        # mesh and function spaces
        self.mesh = Mesh(const.MESH_XML(ref))
        self.V = VectorFunctionSpace(self.mesh, const.V, const.V_DIM)
        self.Q = FunctionSpace(self.mesh, const.Q, const.Q_DIM)
        self.boundaryfunction = MeshFunction("size_t", self.mesh, const.BOUNDARY_XML(self.ref))

        # define attributes v_inner_nodes and v_outer_nodes
        if self.flag == "sim":
            # parts of boundary are added as inner nodes
            self.inner_nodes, self.outer_nodes = inner_outer_nodes(self.V, self.boundaryfunction, \
                                                                   const.ASSEMBLER_COMPRESS_SIM_INNERNODES)
        elif self.flag == "ctrl":
            # GammaRight Boundary Part is added to the inner nodes
            self.inner_nodes, self.outer_nodes = inner_outer_nodes(self.V, self.boundaryfunction, \
                                                                   const.ASSEMBLER_COMPRESS_CTRL_INNERNODES)
        else:
            raise ValueError('unknown flag {s}'.format(self.flag))

        # read matrices and compress system
        names = ["M", "M_BOUNDARY_CTRL", "S", "R", "K", "G", "GT", "B", "C"]
        self.mat = {}
        self.matcps = {}
        for name in names:
            self.mat[name] = scio.mmread(const.ASSEMBLER_NAME_MTX(ref, name, RE))
            self.matcps[name] = compress_mat(self.mat[name], self.inner_nodes)



    def save_mtx(self):
        """save compress system matrices"""

        for name, mat in self.matcps.iteritems():
            if self.flag == "sim":
                filename = self.const.ASSEMBLER_COMPRESS_SIM_NAME_MTX(self.ref, name, self.RE)
            elif self.flag == "ctrl":
                filename = self.const.ASSEMBLER_COMPRESS_CTRL_NAME_MTX(self.ref, name, self.RE)
            else:
                raise ValueError("Unknow flag {0:s}".format(self.flag))

            createdir(filename)
            #with open(filename, "w") as handle:
            #    write_matrix(handle, self.matcps[name], "{0:s},ref={1:d},RE={2:d}, flag={3:s}".format(name, self.ref, self.RE, self.flag))
            #    scio.mmwrite(handle, self.matcps[name])
            if hasattr(mat, "eliminate_zeros"):
                mat.eliminate_zeros()
            write_matrix(filename, mat, "{0:s},ref={1:d},RE={2:d}".format(name, self.ref, self.RE))


        if self.flag == "sim":
            with open(self.const.ASSEMBLER_COMPRESS_SIM_INNERNODES_DAT(self.ref, self.RE), "w") as handle:
                np.savetxt(handle, self.inner_nodes, fmt='%d')

            with open(self.const.ASSEMBLER_COMPRESS_SIM_OUTERNODES_DAT(self.ref, self.RE), "w") as handle:
                np.savetxt(handle, self.outer_nodes, fmt='%d')
        else:
            with open(self.const.ASSEMBLER_COMPRESS_CTRL_INNERNODES_DAT(self.ref, self.RE), "w") as handle:
                np.savetxt(handle, self.inner_nodes, fmt='%d')

            with open(self.const.ASSEMBLER_COMPRESS_CTRL_OUTERNODES_DAT(self.ref, self.RE), "w") as handle:
                np.savetxt(handle, self.outer_nodes, fmt='%d')


