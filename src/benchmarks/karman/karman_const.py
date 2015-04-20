# -*- coding: utf-8 -*-
import dolfin
from dolfin import parameters
from dolfin import Expression
from dolfin import Constant
from dolfin import DirichletBC
import socket
import numpy as np
import os

"""name of scneario"""
NAME = "karman"

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
RECT = {"x0_0": 0, "x1_0": 5.0, "x0_1": 0, "x1_1": 1.0}

"""circular obstacle domain"""
CIRCLE = {"x0": 0.5, "x1": 0.5, "r": 0.10, "fragments": 8}

"""resolution of the macro mesh"""
INITIALRESOLUTION = 12

"""description of upper control part"""
GAMMA_BALL_CTRLUPPER_UPPER_X1 = CIRCLE["x1"] + 6.0/8.0 * CIRCLE["r"]
GAMMA_BALL_CTRLUPPER_LOWER_X1 = CIRCLE["x1"] + 1.0/8.0 * CIRCLE["r"]
GAMMA_BALL_CTRLUPPER_THRESHOLD = 1.2

assert GAMMA_BALL_CTRLUPPER_THRESHOLD > 1, "Constant must be slightly larger than 1"

"""description of lower control part"""
GAMMA_BALL_CTRLLOWER_UPPER_X1 = CIRCLE["x1"] - 1.0/8.0 * CIRCLE["r"]
GAMMA_BALL_CTRLLOWER_LOWER_X1 = CIRCLE["x1"] - 6.0/8.0 * CIRCLE["r"]
GAMMA_BALL_CTRLLOWER_THRESHOLD = 1.2

assert GAMMA_BALL_CTRLLOWER_THRESHOLD > 1, "Constant must be slightly larger than 1"

"""description for boundary obstacle"""
GAMMA_BALL_THRESHOLD = 1.2

assert GAMMA_BALL_THRESHOLD > 1, "Constant must be slightly larger than 1"

"""indices for the boundary parts"""
GAMMA_INNER_INDICES = 0
GAMMA_LEFT_INDICES = 1
GAMMA_LOWER_INDICES = 2
GAMMA_RIGHT_INDICES = 3
GAMMA_UPPER_INDICES = 4
GAMMA_BALL_INDICES = 5
GAMMA_BALL_CTRLLOWER_INDICES = 6
GAMMA_BALL_CTRLUPPER_INDICES = 7

"""Threshold for projection of the boundary of the obstacle"""
GAMMA_BALL_PROJECTION_THRESHOLD = 1.1

assert GAMMA_BALL_PROJECTION_THRESHOLD > 1, "Constant must be slightly larger than 1"


"""define local refinementzone"""
LOCALREFINEMENTS = 1
def LOCALREFINE(p):
    if CIRCLE["x0"] < p.x() < 0.8*(RECT["x1_0"]-RECT["x0_0"]) and \
       0.2*(RECT["x1_1"]-RECT["x0_1"])< p.y() < 0.8*(RECT["x1_1"]-RECT["x0_1"]):
        return True
    return False

def GET_NU(RE):
    """return nu for given RE"""
    # characteristic velocity is 1 (maximum of STATIONARY_UIN)
    # characteristic lenght is 1 (height of the mesh)
    # return float(1.0)/float(RE)
    return Expression("R/RE", R=float(CIRCLE["r"]), RE=float(RE))

def GET_NU_FLOAT(RE):
    return float(CIRCLE["r"])/float(RE)

"""velocity for time stepping"""
U = 1.0

"""Output directory"""
OUTPUTDIR_NAME = "results"

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


"""constants for stationary solver"""
STATIONARY_DIR = "stationary"
STATIONARY_RHS = Constant((0.0, 0.0))
STATIONARY_NEWTON_STEPS = 15
STATIONARY_NEWTON_ABS_TOL = 1e-12
STATIONARY_NEWTON_REL_TOL = 1e-14
STATIONARY_NEWTON_REPORT = False
STATIONARY_UIN = Expression(("1.0/(pow(ye/2.0-ya/2.0,2.0))*(ye-x[1])*(x[1]-ya)", "0.0"), \
                            ya=float(RECT["x0_0"]), ye=float(RECT["x1_1"]))


def STATIONARY_BOUNDARY_CONDITIONS(W, boundaryfunction, const):

    # noslip at boundary parts
    noslip = Constant((0.0, 0.0))

    # define and collect boundary conditions
    bcu = [DirichletBC(W.sub(0), const.STATIONARY_UIN, boundaryfunction, GAMMA_LEFT_INDICES),
           DirichletBC(W.sub(0), noslip, boundaryfunction, GAMMA_LOWER_INDICES),
           DirichletBC(W.sub(0), noslip, boundaryfunction, GAMMA_UPPER_INDICES),
           DirichletBC(W.sub(0), noslip, boundaryfunction, GAMMA_BALL_INDICES),
           DirichletBC(W.sub(0), noslip, boundaryfunction, GAMMA_BALL_CTRLLOWER_INDICES),
           DirichletBC(W.sub(0), noslip, boundaryfunction, GAMMA_BALL_CTRLUPPER_INDICES)]
    return bcu


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
ASSEMBLER_PENALTY_EPS = 1e-6
ASSEMBLER_OBSERVER_POINTS = [(2.5, 0.25), (2.5, 0.75)]

ASSEMBLER_UPPER_CONTROL = Expression(
            ("1.0/pow(phi2/2.0-phi1/2.0, 2.0) * (x[0]-x0)* (atan((x[0]-x0)/(x[1]-y0))-phi1) * (phi2-atan((x[0]-x0)/(x[1]-y0)))",
             "1.0/pow(phi2/2.0-phi1/2.0, 2.0) * (x[1]-y0)* (atan((x[0]-x0)/(x[1]-y0))-phi1) * (phi2-atan((x[0]-x0)/(x[1]-y0)))"),
            phi1=np.pi / 2.0 - np.arccos(1.0 / 8.0),
            phi2=np.pi / 2.0 - np.arccos(6.0 / 8.0),
            x0=CIRCLE["x0"],
            y0=CIRCLE["x1"])

ASSEMBLER_LOWER_CONTROL = Expression(
           ("1.0/pow(phi2/2.0-phi1/2.0, 2.0) * (x[0]-x0)* (atan((x[0]-x0)/(x[1]-y0))-phi1) * (phi2-atan((x[0]-x0)/(x[1]-y0)))",
            "1.0/pow(phi2/2.0-phi1/2.0, 2.0) * (x[1]-y0)* (atan((x[0]-x0)/(x[1]-y0))-phi1) * (phi2-atan((x[0]-x0)/(x[1]-y0)))"),
           phi1=-np.pi / 2.0 + np.arccos(1.0 / 8.0),
           phi2=-np.pi / 2.0 + np.arccos(6.0 / 8.0),
           x0=CIRCLE["x0"],
           y0=CIRCLE["x1"])

ASSEMBLER_BOUNDARY_CONTROLS =[(ASSEMBLER_UPPER_CONTROL, GAMMA_BALL_CTRLUPPER_INDICES, "upper"), \
                              (ASSEMBLER_LOWER_CONTROL, GAMMA_BALL_CTRLLOWER_INDICES, "lower")]

def ASSEMBLER_NAME_MTX(ref, name, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_DIR, parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "{0:s}.mtx.gz".format(name))

def ASSEMBLER_MAT(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_DIR, parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "assembler.mat")


"""constant for compress assembler"""
ASSEMBLER_COMPRESS_SIM_DIR = "assembler_compress_sim"
ASSEMBLER_COMPRESS_SIM_DIRI_ZEROS = [GAMMA_LEFT_INDICES, GAMMA_LOWER_INDICES, \
                                     GAMMA_UPPER_INDICES, GAMMA_BALL_INDICES, \
                                     GAMMA_BALL_CTRLLOWER_INDICES, GAMMA_BALL_CTRLUPPER_INDICES]



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
ASSEMBLER_COMPRESS_CTRL_DIRI_ZEROS = [GAMMA_LEFT_INDICES, GAMMA_LOWER_INDICES, \
                                      GAMMA_UPPER_INDICES, GAMMA_BALL_INDICES]

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
LINEARIZED_SIM_SAVE_PER_S = 5
LINEARIZED_SIM_INFO = 0.05
LINEARIZED_SIM_PERTUBATIONEPS = 0.25
LINEARIZED_SIM_DT = 0.002
LINEARIZED_SIM_STABLE_DT = 1
LINEARIZED_SIM_T = 60
LINEARIZED_SIM_CORRECTION_STEPS = 60    # correction steps for time integration scheme
LINEARIZED_SIM_CORRECTION_RES = 1e-15   # correction residual for time intergraion scheme
#LINEARIZED_SIM_RES = 1e-5               # break if ||u_delta|| smaller this bound
LINEARIZED_SIM_RES = 0               # break if ||u_delta|| smaller this bound
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
BERNOULLI_INSTABLE_RE = 80 # a rough estimation when instabilities occur
BERNOULLI_STRATEGY_1 = {"solver": "scipy", "strategy": "shiftinvert", "sigma": 0.5, "eigenvals": 400, "target": "LM", "tol": 1e-13}
BERNOULLI_STRATEGY_2 = {"solver": "scipy", "strategy": "shiftinvert", "sigma": 1.0, "eigenvals": 4, "target": "LR", "tol": 1e-13}
BERNOULLI_STRATEGY_3 = {"solver": "scipy", "strategy": "moebius", "sigma": 1.0, "tau": -1.0, "eigenvals": 4, "tol": 1e-13,
                        "target": "LM"}

# choose eigenvals in moebius and shiftinvert with LR to that number of instable eigenvalues
# for karman experiments have shown that there are two
# sigma must be positive and tau negative and they must have the same absolute value


def BERNOULLI_FEED0_CPS_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_COMPRESS_CTRL_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "Feed0.mtx.gz")

def BERNOULLI_FEED1_CPS_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_COMPRESS_CTRL_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "Feed1.mtx.gz")

"""constants for lqr solver of navier stokes"""
LQR_DELTA = - 0.02
LQR_NM_OUTPUT = 1
LQR_NM_RES2_SAVE = 1e-1
LQR_NM_RES2 = 5e-10
LQR_NM_REL2_CHANGE = 3e-10
LQR_NM_REL_CHANGE = 3e-10
LQR_NM_MAXIT = 20
LQR_ADI_OUTPUT = 1
LQR_ADI_RES2 = 1e-15
LQR_ADI_MAXIT = 1500
LQR_ADI_REL_CHANGE_TOL = 1e-13
LQR_ADI_ARP_M = 200
LQR_ADI_ARP_P = 200
LQR_ADI_L0 = 50
LQR_PARATYPE = 4 # 4-adaptive V 1 minmax
LQR_NM_GP = 0
LQR_GP = 0
LQR_SAVE_FREQ = 5
LQR_START_CONTROLLING = 0
LQR_INFO = 0.05
if socket.gethostname() == "editha":
    LQR_MEMORY_USAGE = 1
else:
    LQR_MEMORY_USAGE = 0

"""constants for linearized control of navier stokes"""
LINEARIZED_CTRL_DIR = "lqr_ctrl"
LINEARIZED_CTRL_SAVE_PER_S = 5
LINEARIZED_CTRL_INFO = 0.05
LINEARIZED_CTRL_PERTUBATIONEPS = 0.25
LINEARIZED_CTRL_DT = 0.002
LINEARIZED_CTRL_STABLE_DT = 1
LINEARIZED_CTRL_T = 60
LINEARIZED_CTRL_START_CONTROLLING = 0.0
LINEARIZED_CTRL_CORRECTION_STEPS = 60    # correction steps for time integration scheme
LINEARIZED_CTRL_CORRECTION_RES = 1e-15   # correction residual for time intergraion scheme
#LINEARIZED_CTRL_RES = 1e-5              # break if ||u_delta|| smaller this bound
LINEARIZED_CTRL_RES = 0                  # break if ||u_delta|| smaller this bound
LINEARIZED_CTRL_CORRECTION_RES_MOD = 5   # in every 5 steps residual is computed


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
