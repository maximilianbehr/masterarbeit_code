import time
import sys
import os
import numpy as np
import dolfin


class Tee(object):
    def __init__(self, *files):
        self.files = files

    def write(self, obj):
        for f in self.files:
            f.write(obj)

    def flush(self):
        for f in self.files:
            f.flush()


class TeeHandler():
    def __init__(self, file):
        dir = os.path.dirname(file)
        if not os.path.exists(dir):
            os.makedirs(dir)
        self.log = open(file, "w")

    def start(self):
        self.original = sys.stdout
        sys.stdout = Tee(sys.stdout, self.log)

    def stop(self):
        sys.stdout = self.original
        self.log.close()


def deletedir(dir):
    fileList = os.listdir(dir)
    if os.path.exists(dir):
        for fileName in fileList:
            os.remove(os.path.join(dir, fileName))
        os.rmdir(dir)


def gettime():
    return time.strftime("%d/%m/%Y %H:%M:%S")


# CODE von Jan Heiland
# dolfin_navier_scipy
def expand_vp_dolfunc(V=None, Q=None, invinds=None, diribcs=None, vp=None,
                      vc=None, pc=None):
    """expand v [and p] to the dolfin function representation
    Parameters
    ----------
    V : dolfin.VectorFunctionSpace
        FEM space of the velocity
    Q : dolfin.FunctionSpace
        FEM space of the pressure
    invinds : (N,) array
        vector of indices of the velocity nodes
    diribcs : list
        of the (Dirichlet) velocity boundary conditions
    vp : (N+M,1) array, optional
        solution vector of velocity and pressure
    v : (N,1) array, optional
        solution vector of velocity
    p : (M,1) array, optional
        solution vector of pressure
    Returns
    -------
    v : dolfin.Function(V)
        velocity as function
    p : dolfin.Function(Q), optional
        pressure as function
    See Also
    --------
    expand_vecnbc_dolfunc : for a scalar function with multiple bcs
    """

    if vp is not None:
        vc = vp[:len(invinds), :]
        pc = vp[len(invinds):, :]
        p = dolfin.Function(Q)
    elif pc is not None:
        p = dolfin.Function(Q)

    v = dolfin.Function(V)
    ve = np.zeros((V.dim(), 1))

    # fill in the boundary values
    for bc in diribcs:
        bcdict = bc.get_boundary_values()
        ve[bcdict.keys(), 0] = bcdict.values()

    ve[invinds] = vc

    if pc is not None:
        pe = np.vstack([pc, [0]])
        p.vector().set_local(pe)
    else:
        p = None

    v.vector().set_local(ve)

    return v, p