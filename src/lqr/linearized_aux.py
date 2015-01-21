import numpy as np
from dolfin import *


def compress_mat(M, inner_nodes):

    if isinstance(M, np.ndarray):
        if M.shape[0] < M.shape[1]:
            ret = M[:, inner_nodes]
        elif M.shape[0] > M.shape[1]:
            ret = M[inner_nodes, :]
        else:
            raise ValueError('M is dense and fullsize matrix?')
    else:
        if M.shape[0] == M.shape[1]:
            ret = M.tocsr()[inner_nodes, :].tocsc()[:, inner_nodes]
            ret.eliminate_zeros()
        elif M.shape[0] < M.shape[1]:
            ret = M.tocsc()[:, inner_nodes]
            ret.eliminate_zeros()
        else:
            ret = M.tocsr()[inner_nodes, :]
            ret.eliminate_zeros()

    return ret


def inner_outer_nodes(V, innerboundaries):
        """compute return indices of inner and outer nodes,
        function assumes homogenius dirichlet boundary conditions
        on boundary parts not specified by innerboundaries"""

        # collect all boundary parts
        noslip = Constant((0.0, 0.0))
        bcouter = DirichletBC(V, noslip, DomainBoundary())
        outer_nodes = set(bcouter.get_boundary_values().keys())

        # remove boundary parts which should belong to the inner
        for innerboundary in innerboundaries:
            # create abritray dirichlet bc for identifying the nodes
            d = DirichletBC(V, noslip, innerboundary)
            outer_nodes = outer_nodes - set(d.get_boundary_values().keys())


        # inner nodes are all nodes which are no outer nodes
        inner_nodes = set(range(0, V.dim())) - outer_nodes

        return np.array(list(inner_nodes)), np.array(list(outer_nodes))

def u_uncompress(V, u, inner_nodes):
    ret = np.zeros((V.dim(),))
    ret[inner_nodes] = u
    return ret

def u_compress(u, inner_nodes):
    return u[inner_nodes]









