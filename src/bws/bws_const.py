import os
import dolfin
from dolfin import parameters
from dolfin import Constant
from dolfin import DirichletBC
from dolfin import Expression
import socket

DOLFIN_VERSION = dolfin.__version__


"""rectangular domain"""
RECTUPPER = {"x0_0": 0.00, "x1_0": 8.00, "x0_1": 0.30, "x1_1": 1.00}
RECTLOWER = {"x0_0": 2.0, "x1_0": 8.00, "x0_1": 0.00, "x1_1": 0.30}

assert RECTUPPER["x0_0"] < RECTLOWER["x0_0"] < RECTUPPER["x1_0"]
assert RECTLOWER["x0_1"] < RECTUPPER["x0_1"] < RECTUPPER["x1_1"]
assert RECTLOWER["x1_0"] == RECTUPPER["x1_0"]
assert RECTLOWER["x1_1"] == RECTUPPER["x0_1"]

"""resolution of the macro mesh"""
INITIALRESOLUTION = 12

"""indices for the boundary parts"""
GAMMA_INNER_INDICES = 0
GAMMA1_INDICES = 1
GAMMA2_INDICES = 2
GAMMA3_INDICES = 3
GAMMA4_INDICES = 4
GAMMA5_INDICES = 5
GAMMA6_INDICES = 6


"""Output directory"""
OUTPUTDIR_NAME = "results_bws"

def OUTPUTDIR():
    dirname = os.path.join(OUTPUTDIR_NAME, DOLFIN_VERSION)
    host = socket.gethostname()
    if host in ["pc747", "pc633", "pc800"]:
        return os.path.abspath(os.path.join("/scratch/behr/master/", dirname))
    elif host == "jack":
        return os.path.abspath(os.path.join("/home/daniels/PycharmProjects/master", dirname))
    elif host == "editha":
        return os.path.abspath(os.path.join("/scratch/vol1/behr/", dirname))
    elif host == "heinrich.mpi-magdeburg.mpg.de":
        return os.path.abspath(os.path.join("/scratch/behr/", dirname))
    else:
        raise NotImplementedError("Outputpath for {0:s} is not implemented".format(host))


"""Output files for mesh generation"""
def BOUNDARY_PVD(ref):
    return os.path.join(OUTPUTDIR(), "mesh", parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "boundary.pvd")


def BOUNDARY_XDMF(ref):
    return os.path.join(OUTPUTDIR(), "mesh", parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "boundary.xdmf")


def BOUNDARY_XML(ref):
    return os.path.join(OUTPUTDIR(), "mesh", parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "boundary.xml.gz")


def MESH_PVD(ref):
    return os.path.join(OUTPUTDIR(), "mesh", parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "mesh.pvd")


def MESH_XDMF(ref):
    return os.path.join(OUTPUTDIR(), "mesh", parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "mesh.xdmf")


def MESH_XML(ref):
    return os.path.join(OUTPUTDIR(), "mesh", parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "mesh.xml.gz")



"""constants for stationary solver"""
STATIONARY_RHS = Constant((0.0, 0.0))
STATIONARY_V = "CG"
STATIONARY_V_DIM = 2
STATIONARY_Q = "CG"
STATIONARY_Q_DIM = 1
STATIONARY_U0 = Constant((0.0, 0.0))
STATIONARY_UIN = Expression(("1.0/( ((-x0)/2+(x1)/2)*((x1)/2-(x0)/2) )*(x[1]-x0)*(x1-x[1])", "0.0"),
                            x0=RECTUPPER["x0_1"], x1=RECTUPPER["x1_1"])
STATIONARY_NEWTON_STEPS = 40
STATIONARY_NEWTON_ABS_TOL = 1e-12
STATIONARY_NEWTON_REL_TOL = 1e-14


def STATIONARY_BOUNDARY_CONDITIONS(W):
    from src.bws.mesh.bws import Gamma1
    from src.bws.mesh.bws import Gamma2
    from src.bws.mesh.bws import Gamma3
    from src.bws.mesh.bws import Gamma4
    from src.bws.mesh.bws import Gamma5
    from src.bws.mesh.bws import Gamma6


    # inflow profile
    uin = STATIONARY_UIN

    # noslip at boundary parts
    noslip = Constant((0.0, 0.0))

    # define and collect boundary conditions
    bcu = [DirichletBC(W.sub(0), uin, Gamma1()),
           DirichletBC(W.sub(0), noslip, Gamma2()),
           DirichletBC(W.sub(0), noslip, Gamma3()),
           DirichletBC(W.sub(0), noslip, Gamma4()),
           DirichletBC(W.sub(0), noslip, Gamma6())]
    return bcu


def STATIONARY_NU(RE):
    """return nu for given RE"""
    # characteristic velocity is 1 (maximum of STATIONARY_UIN)
    # characteristic lenght is 1 (height of the mesh)
    return Constant(1.0/RE)


def STATIONARY_U_PVD(ref, RE):
    return os.path.join(OUTPUTDIR(), "stationary", parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "u.pvd")


def STATIONARY_U_XML(ref, RE):
    return os.path.join(OUTPUTDIR(), "stationary", parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "u.xml")


def STATIONARY_P_PVD(ref, RE):
    return os.path.join(OUTPUTDIR(), "stationary", parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "p.pvd")


def STATIONARY_P_XML(ref, RE):
    return os.path.join(OUTPUTDIR(), "stationary", parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "p.xml")


def STATIONARY_W_XML(ref, RE):
    return os.path.join(OUTPUTDIR(), "stationary", parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "w.xml")


