import src.karman_const as const
from linearized_aux import compress_mat
from linearized_aux import inner_outer_nodes
from src.mesh.karman import GammaRight
from src.mesh.karman import GammaBallCtrlLower
from src.mesh.karman import GammaBallCtrlUpper
import scipy.io as scio
import numpy as np
import os
from dolfin import *

class CompressAssembler():
    def __init__(self, ref, RE, flag):

        # parameters
        self.ref = ref
        self.RE = RE
        self.flag = flag

        # mesh and function spaces
        self.mesh = Mesh(const.MESH_XML(ref))
        self.V = VectorFunctionSpace(self.mesh, const.LINEARIZED_SIM_V, const.LINEARIZED_SIM_V_DIM)
        self.Q = FunctionSpace(self.mesh, const.LINEARIZED_SIM_Q, const.LINEARIZED_SIM_Q_DIM)

        # define attributes v_inner_nodes and v_outer_nodes
        if self.flag == 'sim':
            # GammaRight Boundary Part is added to the inner nodes
            self.inner_nodes, self.outer_nodes = inner_outer_nodes(self.V, [GammaRight()])
        elif self.flag == 'ctrl':
            # GammaRight Boundary Part is added to the inner nodes
            self.inner_nodes, self.outer_nodes = inner_outer_nodes(self.V, [GammaRight(), GammaBallCtrlLower(), GammaBallCtrlUpper()])
        else:
            raise ValueError('unknown flag {s}'.format(self.flag))

        # read matrices
        self.M = scio.mmread(const.ASSEMBLER_M_MTX(ref, RE))
        self.Mlower = scio.mmread(const.ASSEMBLER_MLOWER_MTX(ref, RE))
        self.Mupper = scio.mmread(const.ASSEMBLER_MUPPER_MTX(ref, RE))
        self.S = scio.mmread(const.ASSEMBLER_S_MTX(ref, RE))
        self.R = scio.mmread(const.ASSEMBLER_R_MTX(ref, RE))
        self.K = scio.mmread(const.ASSEMBLER_K_MTX(ref, RE))
        self.G = scio.mmread(const.ASSEMBLER_G_MTX(ref, RE))
        self.GT = scio.mmread(const.ASSEMBLER_GT_MTX(ref, RE))
        self.B = scio.mmread(const.ASSEMBLER_B_MTX(ref, RE))
        self.C = scio.mmread(const.ASSEMBLER_C_MTX(ref, RE))

        # compress system
        self.Mcps = compress_mat(self.M, self.inner_nodes)
        self.Mcpslower = compress_mat(self.Mlower, self.inner_nodes)
        self.Mcpsupper = compress_mat(self.Mupper, self.inner_nodes)
        self.Scps = compress_mat(self.S, self.inner_nodes)
        self.Rcps = compress_mat(self.R, self.inner_nodes)
        self.Kcps = compress_mat(self.K, self.inner_nodes)
        self.Gcps = compress_mat(self.G, self.inner_nodes)
        self.GTcps = compress_mat(self.GT, self.inner_nodes)
        self.Bcps = compress_mat(self.B, self.inner_nodes)
        self.Ccps = compress_mat(self.C, self.inner_nodes)


    def save(self):
        """save compress system matrices"""

        matrices = [self.Mcps,
                    self.Mcpslower,
                    self.Mcpsupper,
                    self.Scps,
                    self.Rcps,
                    self.Kcps,
                    self.Gcps,
                    self.GTcps,
                    self.Bcps,
                    self.Ccps]

        if self.flag == 'sim':
            files = [const.ASSEMBLER_COMPRESS_SIM_M_MTX,
                     const.ASSEMBLER_COMPRESS_SIM_MLOWER_MTX,
                     const.ASSEMBLER_COMPRESS_SIM_MUPPER_MTX,
                     const.ASSEMBLER_COMPRESS_SIM_S_MTX,
                     const.ASSEMBLER_COMPRESS_SIM_R_MTX,
                     const.ASSEMBLER_COMPRESS_SIM_K_MTX,
                     const.ASSEMBLER_COMPRESS_SIM_G_MTX,
                     const.ASSEMBLER_COMPRESS_SIM_GT_MTX,
                     const.ASSEMBLER_COMPRESS_SIM_B_MTX,
                     const.ASSEMBLER_COMPRESS_SIM_C_MTX]
        elif self.flag == 'ctrl':
            files = [const.ASSEMBLER_COMPRESS_CTRL_M_MTX,
                     const.ASSEMBLER_COMPRESS_CTRL_MLOWER_MTX,
                     const.ASSEMBLER_COMPRESS_CTRL_MUPPER_MTX,
                     const.ASSEMBLER_COMPRESS_CTRL_S_MTX,
                     const.ASSEMBLER_COMPRESS_CTRL_R_MTX,
                     const.ASSEMBLER_COMPRESS_CTRL_K_MTX,
                     const.ASSEMBLER_COMPRESS_CTRL_G_MTX,
                     const.ASSEMBLER_COMPRESS_CTRL_GT_MTX,
                     const.ASSEMBLER_COMPRESS_CTRL_B_MTX,
                     const.ASSEMBLER_COMPRESS_CTRL_C_MTX]

        if not os.path.exists(os.path.dirname(files[0](self.ref, self.RE))):
            os.makedirs(os.path.dirname(files[0](self.ref, self.RE)))

        for idx in range(len(matrices)):
            with open(files[idx](self.ref, self.RE), "w") as handle:
                scio.mmwrite(handle, matrices[idx])

        if self.flag == "sim":
            with open(const.ASSEMBLER_COMPRESS_SIM_INNER_NODES(self.ref, self.RE), "w") as handle:
                np.savetxt(handle, self.inner_nodes, fmt='%d')

            with open(const.ASSEMBLER_COMPRESS_SIM_OUTER_NODES(self.ref, self.RE), "w") as handle:
                np.savetxt(handle, self.outer_nodes, fmt='%d')
        else:
            with open(const.ASSEMBLER_COMPRESS_CTRL_INNER_NODES(self.ref, self.RE), "w") as handle:
                np.savetxt(handle, self.inner_nodes, fmt='%d')

            with open(const.ASSEMBLER_COMPRESS_CTRL_OUTER_NODES(self.ref, self.RE), "w") as handle:
                np.savetxt(handle, self.outer_nodes, fmt='%d')


