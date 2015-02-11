import os
import dolfin
from dolfin import parameters
from dolfin import Constant
from dolfin import DirichletBC
from dolfin import Expression
import socket
import numpy as np

DOLFIN_VERSION = dolfin.__version__


"""rectangular domain"""
# parameter for adapting size of backward facing step
MODELHEIGHT = 1.0
RECTLOWER = {"x0_0": 5*MODELHEIGHT, "x0_1": 0.00, "x1_0": 25*MODELHEIGHT, "x1_1": MODELHEIGHT}
RECTUPPER = {"x0_0": 0.00, "x0_1": MODELHEIGHT, "x1_0": 25*MODELHEIGHT, "x1_1": 5*MODELHEIGHT}
RECTROTATE = {"diag": 0.25*MODELHEIGHT, "angle": np.pi/4}

assert RECTUPPER["x0_0"] < RECTLOWER["x0_0"] < RECTUPPER["x1_0"]
assert RECTLOWER["x0_1"] < RECTUPPER["x0_1"] < RECTUPPER["x1_1"]
assert RECTLOWER["x1_0"] == RECTUPPER["x1_0"]
assert RECTLOWER["x1_1"] == RECTUPPER["x0_1"]

"""define local refinementzone"""
LOCALREFINEMENTS = 2
def LOCALREFINE(p):
    if 3*MODELHEIGHT < p.x() < 20*MODELHEIGHT and p.y() < 1.5*MODELHEIGHT:
        return True
    return False

"""resolution of the macro mesh"""
INITIALRESOLUTION = 12

"""indices for the boundary parts"""
CONTROLRADIUS = 0.25*MODELHEIGHT
GAMMA_INNER_INDICES = 0
GAMMA1_INDICES = 1
GAMMA2_INDICES = 2
GAMMA3_INDICES = 3
GAMMA4_INDICES = 4
GAMMA5_INDICES = 5
GAMMA6_INDICES = 6
GAMMA7_INDICES = 7
GAMMA8_INDICES = 8
GAMMASHEAR1_INDICES = 9
GAMMASHEAR2_INDICES = 10
GAMMASHEAR3_INDICES = 11


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


def SHEAR_PVD(ref, num):
    return os.path.join(OUTPUTDIR(), "mesh", parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "shear{0:d}.pvd".format(num))


def SHEAR_XDMF(ref, num):
    return os.path.join(OUTPUTDIR(), "mesh", parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "shear{0:d}.xdmf".format(num))


def SHEAR_XML(ref, num):
    return os.path.join(OUTPUTDIR(), "mesh", parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "shear{0:d}.xml.gz".format(num))


"""constants for stationary solver"""
STATIONARY_RHS = Constant((0.0, 0.0))
STATIONARY_V = "CG"
STATIONARY_V_DIM = 2
STATIONARY_Q = "CG"
STATIONARY_Q_DIM = 1
STATIONARY_U0 = Constant((0.0, 0.0))
STATIONARY_UIN = Expression(("1.0/pow((x1)/2.0-(x0)/2.0,2) *(x[1]-x0)*(x1-x[1])", "0.0"),
                            x0=RECTUPPER["x0_1"], x1=RECTUPPER["x1_1"])
STATIONARY_NEWTON_STEPS = 15
STATIONARY_NEWTON_ABS_TOL = 1e-12
STATIONARY_NEWTON_REL_TOL = 1e-14

l = RECTROTATE["diag"]/np.sqrt(2)
# middlepoint of schraege
x_m = RECTLOWER["x0_0"] - l*np.cos(RECTROTATE["angle"])
y_m = RECTUPPER["x0_1"] - l*np.sin(RECTROTATE["angle"])
# anfangspoints of schraege
x_a = x_m - CONTROLRADIUS/2.0*(np.sin(RECTROTATE["angle"]))
y_a = y_m + CONTROLRADIUS/2.0*(np.cos(RECTROTATE["angle"]))
# endpoint of schraege
x_e = x_m + CONTROLRADIUS/2.0*(np.sin(RECTROTATE["angle"]))
y_e = y_m - CONTROLRADIUS/2.0*(np.cos(RECTROTATE["angle"]))

STATIONARY_CONTROL_LAMBDA = 0.2
# print "pm=[{0:e},{1:e}];".format(x_m, y_m)
# print "pa=[{0:e},{1:e}];".format(x_a, y_a)
# print "pe=[{0:e},{1:e}];".format(x_e, y_e)

STATIONARY_CONTROL_UPPER = Expression(("lam*sqrt(pow(x[0]-xa,2)+pow(x[1]-ya,2))*sqrt(pow(x[0]-xe,2)+pow(x[1]-ye,2))/pow(r/2.0,2.0)*cos(phi)", \
                                       "lam*sqrt(pow(x[0]-xa,2)+pow(x[1]-ya,2))*sqrt(pow(x[0]-xe,2)+pow(x[1]-ye,2))/pow(r/2.0,2.0)*sin(phi)"), \
                                       r=CONTROLRADIUS, phi=RECTROTATE["angle"], xa=x_a, ya=y_a, xe=x_e, ye=y_e, \
                                       lam=STATIONARY_CONTROL_LAMBDA)

#STATIONARY_CONTROL_LOWER = Expression(("lam*(-1.0)*(x[1]*(r-x[1]))/pow(r/2.0,2.0)", "0.0"), r=CONTROLRADIUS, lam=STATIONARY_CONTROL_LAMBDA)
STATIONARY_CONTROL_LOWER = Expression(("-2*(-1.0)*(x[1]*(r-x[1]))/pow(r/2.0,2.0)", "0.0"), r=CONTROLRADIUS, lam=STATIONARY_CONTROL_LAMBDA)


def STATIONARY_BOUNDARY_CONDITIONS(W, boundaryfunction):
    from src.bws.mesh.bws import Gamma1
    from src.bws.mesh.bws import Gamma2
    from src.bws.mesh.bws import Gamma3
    from src.bws.mesh.bws import Gamma4
    from src.bws.mesh.bws import Gamma5
    from src.bws.mesh.bws import Gamma6
    from src.bws.mesh.bws import Gamma7
    from src.bws.mesh.bws import Gamma8


    # inflow profile
    uin = STATIONARY_UIN

    # control upper in
    ctrl_upper = STATIONARY_CONTROL_UPPER

    # control lower in
    ctrl_lower = STATIONARY_CONTROL_LOWER


    # noslip at boundary parts
    noslip = Constant((0.0, 0.0))

    # define and collect boundary conditions
    bcu = [DirichletBC(W.sub(0), uin, boundaryfunction, Gamma1().index),
           DirichletBC(W.sub(0), noslip, boundaryfunction, Gamma2().index),

           DirichletBC(W.sub(0), ctrl_upper, boundaryfunction, Gamma3().index),
           # DirichletBC(W.sub(0), noslip, boundaryfunction, Gamma3().index),

           DirichletBC(W.sub(0), noslip, boundaryfunction, Gamma4().index),

           DirichletBC(W.sub(0), ctrl_lower, boundaryfunction, Gamma5().index),
           #DirichletBC(W.sub(0), uin, boundaryfunction, Gamma5().index),

           DirichletBC(W.sub(0), noslip, boundaryfunction, Gamma6().index),
           DirichletBC(W.sub(0), noslip, boundaryfunction, Gamma8().index)]
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


