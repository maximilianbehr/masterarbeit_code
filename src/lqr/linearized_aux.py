import numpy as np
from dolfin import *
import scipy.linalg  as scla

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


def inner_outer_nodes(V, boundaryfunction, diri_zeros):
        """compute return indices of inner and outer nodes,
        function assumes homogenius dirichlet boundary conditions
        on boundary parts not specified by innerboundaries"""

        # collect all boundary parts
        noslip = Constant((0.0, 0.0))

        # remove boundary parts which should belong to the inner
        diri_nodes =[]
        for diri_zero in diri_zeros:
            # create abritray dirichlet bc for identifying the nodes
            d = DirichletBC(V, noslip, boundaryfunction, diri_zero)
            diri_nodes.extend(d.get_boundary_values().keys())

        diri_nodes = set(diri_nodes)
        # diri nodes are all nodes which are no outer nodes
        inner_nodes = set(range(0, V.dim())) - diri_nodes

        diri_nodes = np.array(list(diri_nodes))
        diri_nodes.sort()
        inner_nodes = np.array(list(inner_nodes))
        inner_nodes.sort()
        return diri_nodes, inner_nodes

def smw(smvlu_piv,Msolver,B,K,rhs):
    "perform inv(M-B*K.T)*rhs"
    Minvrhs = Msolver.solve(rhs)
    temp = np.dot(K.T, Minvrhs)
    # Minvb = Msolver.solve(B)
    # temp2 = np.dot(K.T, Minvb)
    # temp2 = np.eye(temp2.shape[0]) - temp2
    temp = scla.lu_solve(smvlu_piv, temp)
    # temp = np.linalg.solve(temp2, temp)
    temp = np.dot(B, temp)
    return Minvrhs + Msolver.solve(temp)


def stable_timestep(T, nu, U, h):
        "Return time step."

        # variant 1
        # dt =  0.25*h**2 / (U*(nu + h*U))

        # variant 2
        dt = 0.2*(float(h) / float(U))
        n = int(float(T) / float(dt) + 1.0)
        dt = float(T) / float(n)
        print "Computing time step according to stability criteria"
        return dt










