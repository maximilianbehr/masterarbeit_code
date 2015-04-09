import os
import dolfin
from dolfin import parameters
from dolfin import Constant
from dolfin import DirichletBC
from dolfin import Expression
import socket
import numpy as np

"""name of scenario"""
NAME = "bws"

"""dolfin version"""
DOLFIN_VERSION = dolfin.__version__

"""label and name for pvd files"""
PVD_U_LABEL_NAME = ("v", "velocity")
PVD_P_LABEL_NAME = ("p", "pressure")
PVD_W_LABEL_NAME = ("v x p", "velocity x pressure")
PVD_MESH_LABEL_NAME = ("mesh", "mesh")

"""function spaces"""
V = "CG"
V_DIM = 2
Q = "CG"
Q_DIM = 1

"""turn dof reordering off"""
parameters["reorder_dofs_serial"] = False

"""set refinement algorithm"""
parameters["refinement_algorithm"] = "regular_cut"

"""Optimzation parameters"""
parameters["form_compiler"]["optimize"] = True
parameters["form_compiler"]["cpp_optimize"] = True
parameters["krylov_solver"]["absolute_tolerance"] = 1e-25
parameters["krylov_solver"]["relative_tolerance"] = 1e-14
parameters["krylov_solver"]["monitor_convergence"] = False


"""rectangular domain"""
# parameter for adapting size of backward facing step
MODELHEIGHT = 1.0
RECTLOWER = {"x0_0": 5.0*MODELHEIGHT, "x0_1": 0.00, "x1_0": 25*MODELHEIGHT, "x1_1": MODELHEIGHT}
RECTUPPER = {"x0_0": 0.00, "x0_1": 1*MODELHEIGHT, "x1_0": 25*MODELHEIGHT, "x1_1": 2.0*MODELHEIGHT}

assert RECTUPPER["x0_0"] < RECTLOWER["x0_0"] < RECTUPPER["x1_0"]
assert RECTLOWER["x0_1"] < RECTUPPER["x0_1"] < RECTUPPER["x1_1"]
assert RECTLOWER["x1_0"] == RECTUPPER["x1_0"]
assert RECTLOWER["x1_1"] == RECTUPPER["x0_1"]

"""define local refinementzone"""
LOCALREFINEMENTS = 1
def LOCALREFINE(p):
    if ((RECTUPPER["x1_0"]*0.2) < p.x() < 0.4*RECTUPPER["x1_0"]) and p.y() <= 1.0*RECTUPPER["x1_1"]:
        return True
    return False

"""resolution of the macro mesh"""
INITIALRESOLUTION = 80

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
    if host in ["pc747", "pc633", "pc800", "pc731", "adelheid", "heinrich.mpi-magdeburg.mpg.de"]:
        return os.path.abspath(os.path.join("/scratch/behr/master/", dirname))
    elif host == "jack":
        return os.path.abspath(os.path.join("/home/daniels/PycharmProjects/master", dirname))
    elif host == "editha":
        return os.path.abspath(os.path.join("/scratch/vol1/behr/", dirname))
    else:
        raise NotImplementedError("Outputpath for {0:s} is not implemented".format(host))


"""Output files for mesh generation"""
def BOUNDARY_PVD(ref):
    return os.path.join(OUTPUTDIR(), "mesh", parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "boundary.pvd")


def BOUNDARY_XML(ref):
    return os.path.join(OUTPUTDIR(), "mesh", parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "boundary.xml.gz")


def MESH_PVD(ref):
    return os.path.join(OUTPUTDIR(), "mesh", parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "mesh.pvd")


def MESH_XML(ref):
    return os.path.join(OUTPUTDIR(), "mesh", parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "mesh.xml.gz")


def SHEAR_PVD(ref, num):
    return os.path.join(OUTPUTDIR(), "mesh", parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "shear{0:d}.pvd".format(num))


def SHEAR_XML(ref, num):
    return os.path.join(OUTPUTDIR(), "mesh", parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "shear{0:d}.xml.gz".format(num))


"""constants for stationary solver"""
STATIONARY_DIR = "stationary"
STATIONARY_RHS = Constant((0.0, 0.0))
STATIONARY_NEWTON_STEPS = 15
STATIONARY_NEWTON_ABS_TOL = 1e-12
STATIONARY_NEWTON_REL_TOL = 1e-14
STATIONARY_NEWTON_REPORT = False
STATIONARY_LAM = 0.0
STATIONARY_CONTROL_UPPER = Expression(("lam*( 1.0)*((h-x[1])*(x[1]-(h-r)))/pow(r/2.0,2.0)", "0.0"), r=CONTROLRADIUS, h=RECTUPPER["x0_1"], lam=STATIONARY_LAM)
STATIONARY_CONTROL_LOWER = Expression(("lam*(-1.0)*(x[1]*(r-x[1]))/pow(r/2.0,2.0)", "0.0"), r=CONTROLRADIUS, lam=STATIONARY_LAM)
STATIONARY_UIN = Expression(("1.0/pow((x1)/2.0-(x0)/2.0,2) *(x[1]-x0)*(x1-x[1])", "0.0"), x0=RECTUPPER["x0_1"], x1=RECTUPPER["x1_1"])



def STATIONARY_BOUNDARY_CONDITIONS(W, boundaryfunction, const):

    # noslip at boundary parts
    noslip = Constant((0.0, 0.0))

    # update boundary conditions
    const.STATIONARY_CONTROL_UPPER.lam = const.STATIONARY_LAM
    const.STATIONARY_CONTROL_LOWER.lam = const.STATIONARY_LAM

    # define and collect boundary conditions
    bcu = [DirichletBC(W.sub(0), STATIONARY_UIN, boundaryfunction, GAMMA1_INDICES),
           DirichletBC(W.sub(0), noslip, boundaryfunction, GAMMA2_INDICES),
           DirichletBC(W.sub(0), const.STATIONARY_CONTROL_UPPER, boundaryfunction, GAMMA3_INDICES),
           DirichletBC(W.sub(0), noslip, boundaryfunction, GAMMA4_INDICES),
           DirichletBC(W.sub(0), const.STATIONARY_CONTROL_LOWER, boundaryfunction, GAMMA5_INDICES),
           DirichletBC(W.sub(0), noslip, boundaryfunction, GAMMA6_INDICES),
           DirichletBC(W.sub(0), noslip, boundaryfunction, GAMMA8_INDICES)]
    return bcu

def GET_NU(RE):
    """return nu for given RE"""
    # characteristic velocity is 1 (maximum of STATIONARY_UIN)
    # characteristic lenght is 1 (height step)
    return Expression("M/RE", M=float(RECTLOWER["x1_1"]-RECTLOWER["x0_1"]), RE=float(RE))

def GET_NU_FLOAT(RE):
    """return nu for given RE"""
    return float(RECTLOWER["x1_1"]-RECTLOWER["x0_1"])/float(RE)

"""velocity for time stepping"""
U = 1.0

def STATIONARY_U_PVD(ref, RE):
    return os.path.join(OUTPUTDIR(), STATIONARY_DIR, parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "u.pvd")


def STATIONARY_U_XML(ref, RE):
    return os.path.join(OUTPUTDIR(), STATIONARY_DIR, parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "u.xml.gz")


def STATIONARY_P_PVD(ref, RE):
    return os.path.join(OUTPUTDIR(), STATIONARY_DIR, parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "p.pvd")


def STATIONARY_P_XML(ref, RE):
    return os.path.join(OUTPUTDIR(), STATIONARY_DIR, parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "p.xml.gz")


def STATIONARY_W_XML(ref, RE):
    return os.path.join(OUTPUTDIR(), STATIONARY_DIR, parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "w.xml.gz")


"""constant for assembler"""
ASSEMBLER_DIR = "assembler"
ASSEMBLER_PENALTY_EPS = 0.001
ASSEMBLER_OBSERVER_POINTS = [(8*MODELHEIGHT, 0.5*MODELHEIGHT), (8*MODELHEIGHT, 1.5*MODELHEIGHT)]


e1 = Expression(("( 1.0)*((h-x[1])*(x[1]-(h-r)))/pow(r/2.0,2.0)", "0.0"), r=CONTROLRADIUS, h=RECTUPPER["x0_1"])
e2 = Expression(("(-1.0)*(x[1]*(r-x[1]))/pow(r/2.0,2.0)", "0.0"), r=CONTROLRADIUS)
ASSEMBLER_BOUNDARY_CONTROLS =[(e1, GAMMA3_INDICES, "upper"), (e2, GAMMA5_INDICES, "lower")]

def ASSEMBLER_NAME_MTX(ref, name, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_DIR, parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "{0:s}.mtx.gz".format(name))

def ASSEMBLER_MAT(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_DIR, parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "assembler.mat")


"""constant for compress assembler simulation"""
ASSEMBLER_COMPRESS_SIM_DIR = "assembler_compress_sim"
ASSEMBLER_COMPRESS_SIM_INNERNODES = [GAMMA7_INDICES]
ASSEMBLER_COMPRESS_SIM_DIRI_ZEROS = [GAMMA1_INDICES, GAMMA2_INDICES, GAMMA3_INDICES,\
                                     GAMMA4_INDICES, GAMMA5_INDICES, GAMMA6_INDICES,\
                                     GAMMA8_INDICES]

def ASSEMBLER_COMPRESS_SIM_NAME_MTX(ref, name, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_COMPRESS_SIM_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "{0:s}.mtx.gz".format(name))

def ASSEMBLER_COMPRESS_SIM_INNERNODES_DAT(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_COMPRESS_SIM_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "innernodes.dat")

def ASSEMBLER_COMPRESS_SIM_OUTERNODES_DAT(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_COMPRESS_SIM_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "outernodes.dat")

"""constant for compress assembler control"""
ASSEMBLER_COMPRESS_CTRL_DIR = "assembler_compress_ctrl"
ASSEMBLER_COMPRESS_CTRL_DIRI_ZEROS = [GAMMA1_INDICES, GAMMA2_INDICES, GAMMA4_INDICES, \
                                      GAMMA6_INDICES, GAMMA8_INDICES]

def ASSEMBLER_COMPRESS_CTRL_NAME_MTX(ref, name, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_COMPRESS_CTRL_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "{0:s}.mtx.gz".format(name))

def ASSEMBLER_COMPRESS_CTRL_INNERNODES_DAT(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_COMPRESS_CTRL_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "innernodes.dat")

def ASSEMBLER_COMPRESS_CTRL_OUTERNODES_DAT(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_COMPRESS_CTRL_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "outernodes.dat")


"""constant for simulation of linearized navier stokes"""
LINEARIZED_SIM_DIR = "lqr_sim"
LINEARIZED_SIM_SAVE_PER_S = 10           # pictures per second
LINEARIZED_SIM_INFO = 0.05
LINEARIZED_SIM_PERTUBATIONEPS = 0.1
LINEARIZED_SIM_DT = 0.01
LINEARIZED_SIM_T = 90
#LINEARIZED_SIM_RES = 1e-5               # break if ||u_delta|| smaller this bound
LINEARIZED_SIM_RES = 0                  # break if ||u_delta|| smaller this bound
LINEARIZED_SIM_STABLE_DT = 1
LINEARIZED_SIM_CORRECTION_STEPS = 60    # correction steps for time integration scheme
LINEARIZED_SIM_CORRECTION_RES = 1e-15   # correction residual for time intergraion scheme
LINEARIZED_SIM_CORRECTION_RES_MOD = 5   # in every 5 steps residual is computed


def LINEARIZED_SIM_U_PVD(ref, RE):
    return os.path.join(OUTPUTDIR(), LINEARIZED_SIM_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "u.pvd")

def LINEARIZED_SIM_U_DELTA_PVD(ref, RE):
    return os.path.join(OUTPUTDIR(), LINEARIZED_SIM_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "u_delta.pvd")

def LINEARIZED_SIM_LOG(ref, RE):
    return os.path.join(OUTPUTDIR(), LINEARIZED_SIM_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "log.txt")

"""constants for bernoulli"""
BERNOULLI_MAXIT = 50
BERNOULLI_STRATEGY = {"strategy": "shiftinvert", "sigma": 1.0, "eigenvals": 1000, "delta": -0.02}
BERNOULLI_STRATEGY = {"strategy": "moebius", "sigma": 1.0, "tau": 1.0, "eigenvals": 1000, "delta": -0.02}


def BERNOULLI_FEED0_CPS_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_COMPRESS_CTRL_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "Feed0.mtx.gz")

def BERNOULLI_FEED1_CPS_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_COMPRESS_CTRL_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "Feed1.mtx.gz")

"""constants for lqr solver of navier stokes"""
LQR_DELTA = - 0.02
LQR_NM_OUTPUT = 1
LQR_NM_RES2_SAVE = 1e-5
LQR_NM_RES2 = 5e-10
LQR_NM_REL2_CHANGE = 3e-10
LQR_NM_REL_CHANGE = 3e-10
LQR_NM_MAXIT = 20
LQR_ADI_OUTPUT = 0
LQR_ADI_RES2 = 1e-16
LQR_ADI_MAXIT = 8000
LQR_ADI_REL_CHANGE_TOL = 1e-13
LQR_ADI_ARP_M = 40
LQR_ADI_ARP_P = 40
LQR_ADI_L0 = 25
LQR_SAVE_FREQ = 5
LQR_START_CONTROLLING = 0
LQR_INFO = 0.05

if socket.gethostname() == "editha":
    LQR_MEMORY_USAGE = 1
else:
    LQR_MEMORY_USAGE = 0

"""constants for linearized control of navier stokes"""
LINEARIZED_CTRL_DIR = "lqr_ctrl"
LINEARIZED_CTRL_SAVE_PER_S = 10
LINEARIZED_CTRL_INFO = 0.05
LINEARIZED_CTRL_PERTUBATIONEPS = 0.25
LINEARIZED_CTRL_DT = 0.01
LINEARIZED_CTRL_T = 90
LINEARIZED_CTRL_STABLE_DT = 1
LINEARIZED_CTRL_START_CONTROLLING = 0.0
LINEARIZED_CTRL_CORRECTION_STEPS = 90       # correction steps for time integration scheme
LINEARIZED_CTRL_CORRECTION_RES = 1e-15      # correction residual for time intergraion scheme
#LINEARIZED_CTRL_RES = 1e-5                  # break if ||u_delta|| smaller this bound
LINEARIZED_CTRL_RES = 0                  # break if ||u_delta|| smaller this bound
LINEARIZED_CTRL_CORRECTION_RES_MOD = 5      # in every 5 steps residual is computed


def LINEARIZED_CTRL_U_PVD(ref, RE):
    return os.path.join(OUTPUTDIR(), LINEARIZED_CTRL_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "u.pvd")

def LINEARIZED_CTRL_U_DELTA_PVD(ref, RE):
    return os.path.join(OUTPUTDIR(), LINEARIZED_CTRL_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "u_delta.pvd")

def LINEARIZED_CTRL_LOG(ref, RE):
    return os.path.join(OUTPUTDIR(), LINEARIZED_CTRL_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "log.txt")


"""constant for eigenvalues"""
EIGEN_DIR = "eigen"

def EIGEN_SYS_CPS_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), EIGEN_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "eig_sys.mtx.gz")

def EIGEN_SYS_CPS_PLOT(ref, RE, ending):
    return os.path.join(OUTPUTDIR(), EIGEN_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "eig_sys.{0:s}".format(ending))

def EIGEN_RIC_CPS_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), EIGEN_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "eig_ric.mtx.gz")

def EIGEN_RIC_CPS_PLOT(ref, RE, ending):
    return os.path.join(OUTPUTDIR(), EIGEN_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "eig_ric.{0:s}".format(ending))

def EIGEN_BER_CPS_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), EIGEN_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "eig_ber.mtx.gz")

def EIGEN_BER_CPS_PLOT(ref, RE, ending):
    return os.path.join(OUTPUTDIR(), EIGEN_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "eig_ber.{0:s}".format(ending))


"""constant for plotter"""
PLOTTER_DIR = "plotter"
def PLOTTER_LINEARIZED_SIM_LOG(ref, num, ending):
    return os.path.join(OUTPUTDIR(), PLOTTER_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "sim{0:d}.".format(num)+ending)

def PLOTTER_LINEARIZED_CTRL_LOG(ref, num, ending):
    return os.path.join(OUTPUTDIR(), PLOTTER_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "ctrl{0:d}.".format(num)+ending)

def PLOTTER_LINEARIZED_CTRL_CONTROLS_LOG(ref, num, ending):
    return os.path.join(OUTPUTDIR(), PLOTTER_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "ctrl_controls{0:d}.".format(num)+ending)

def PLOTTER_LINEARIZED_CTRL_OUTPUT_LOG(ref, num, ending):
    return os.path.join(OUTPUTDIR(), PLOTTER_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "ctrl_outputs{0:d}.".format(num)+ending)

