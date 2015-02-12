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

"""turn dof reordering off"""
parameters["reorder_dofs_serial"] = False

"""set refinement algorithm"""
#parameters["refinement_algorithm"] = "recursive_bisection"
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
CIRCLE = {"x0": 0.5, "x1": 0.5, "r": 0.09, "fragments": 8}

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

"""Output directory"""
OUTPUTDIR_NAME = "results"

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
STATIONARY_UIN = Expression(("4*(1-x[1])*x[1]", "0.0"))
STATIONARY_NEWTON_STEPS = 40
STATIONARY_NEWTON_ABS_TOL = 1e-12
STATIONARY_NEWTON_REL_TOL = 1e-14


def STATIONARY_BOUNDARY_CONDITIONS(W, boundaryfunction):
    from src.benchmarks.karman.mesh.karman import GammaLeft
    from src.benchmarks.karman.mesh.karman import GammaLower
    from src.benchmarks.karman.mesh.karman import GammaUpper
    from src.benchmarks.karman.mesh.karman import GammaBall
    from src.benchmarks.karman.mesh.karman import GammaBallCtrlLower
    from src.benchmarks.karman.mesh.karman import GammaBallCtrlUpper

    # inflow profile
    uin = STATIONARY_UIN

    # noslip at boundary parts
    noslip = Constant((0.0, 0.0))

    # define and collect boundary conditions
    bcu = [DirichletBC(W.sub(0), uin, boundaryfunction, GammaLeft().index),
           DirichletBC(W.sub(0), noslip, boundaryfunction, GammaLower().index),
           DirichletBC(W.sub(0), noslip, boundaryfunction, GammaUpper().index),
           DirichletBC(W.sub(0), noslip, boundaryfunction, GammaBall().index),
           DirichletBC(W.sub(0), noslip, boundaryfunction, GammaBallCtrlLower().index),
           DirichletBC(W.sub(0), noslip, boundaryfunction, GammaBallCtrlUpper().index)]
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

"""constants for instationary solvers"""
INSTATIONARY_RHS = Constant((0.0, 0.0))
INSTATIONARY_V = "CG"
INSTATIONARY_V_DIM = 2
INSTATIONARY_Q = "CG"
INSTATIONARY_Q_DIM = 2
INSTATIONARY_U0 = Constant((0.0, 0.0))
INSTATIONARY_P0 = Constant(0.0)
INSTATIONARY_UIN_MAX = 4
INSTATIONARY_SAVE_FREQ = 5

def INSTATIONARY_UIN(V,Q,t):
    return Expression(("4*(1-x[1])*x[1]*(t/(1.0+t))", "0.0"), t=t)


def INSTATIONARY_NU(RE):
    """return nu for given RE"""
    # characteristic velocity is 1 (maximum of STATIONARY_UIN)
    # characteristic lenght is 1 (height of the mesh)
    return Constant(1.0/RE)


def INSTATIONARY_U_PVD(ref, RE, solver):
    return os.path.join(OUTPUTDIR(), solver, parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "u.pvd")


def INSTATIONARY_P_PVD(ref, RE, solver):
    return os.path.join(OUTPUTDIR(), solver, parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "p.pvd")


"""constant for assembler"""
ASSEMBLER_DIR = "assembler"
ASSEMBLER_V = "CG"
ASSEMBLER_V_DIM = 2
ASSEMBLER_Q = "CG"
ASSEMBLER_Q_DIM = 1
ASSEMBLER_PENALTY_EPS = 0.001
ASSEMBLER_OBSERVER_POINT1_X = 2.5
ASSEMBLER_OBSERVER_POINT1_Y = 0.25
ASSEMBLER_OBSERVER_POINT2_X = 2.5
ASSEMBLER_OBSERVER_POINT2_Y = 0.75
#ASSEMBLER_UPPER_CONTROL = Expression(("1/r * (x[0]-x0)", "1/r * (x[1]-y0)"), r=CIRCLE["r"], x0=CIRCLE["x0"], y0=CIRCLE["x1"])
#ASSEMBLER_LOWER_CONTROL = Expression(("1/r * (x[0]-x0)", "1/r * (x[1]-y0)"), r=CIRCLE["r"], x0=CIRCLE["x0"], y0=CIRCLE["x1"])

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

def ASSEMBLER_M_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_DIR, parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "M.mtx")


def ASSEMBLER_MLOWER_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_DIR, parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "Mlower.mtx")


def ASSEMBLER_MUPPER_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_DIR, parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "Mupper.mtx")


def ASSEMBLER_K_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_DIR, parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "K.mtx")


def ASSEMBLER_R_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_DIR, parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "R.mtx")


def ASSEMBLER_S_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_DIR, parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "S.mtx")


def ASSEMBLER_G_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_DIR, parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "G.mtx")


def ASSEMBLER_GT_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_DIR, parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "GT.mtx")


def ASSEMBLER_BLOWER_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_DIR, parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "Blower.mtx")


def ASSEMBLER_BUPPER_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_DIR, parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "Bupper.mtx")


def ASSEMBLER_B_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_DIR, parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "B.mtx")


def ASSEMBLER_C_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_DIR, parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "C.mtx")


def ASSEMBLER_MAT(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_DIR, parameters["refinement_algorithm"], "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "assembler.mtx")

"""constant for compress assembler"""
ASSEMBLER_COMPRESS_SIM_DIR = "assembler_compress_sim"

def ASSEMBLER_COMPRESS_SIM_M_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_COMPRESS_SIM_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "Mcps.mtx")


def ASSEMBLER_COMPRESS_SIM_MLOWER_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_COMPRESS_SIM_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "Mlowercps.mtx")


def ASSEMBLER_COMPRESS_SIM_MUPPER_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_COMPRESS_SIM_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "Muppercps.mtx")


def ASSEMBLER_COMPRESS_SIM_S_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_COMPRESS_SIM_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE),  "Scps.mtx")


def ASSEMBLER_COMPRESS_SIM_R_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_COMPRESS_SIM_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "Rcps.mtx")


def ASSEMBLER_COMPRESS_SIM_K_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_COMPRESS_SIM_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "Kcps.mtx")


def ASSEMBLER_COMPRESS_SIM_G_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_COMPRESS_SIM_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "Gcps.mtx")


def ASSEMBLER_COMPRESS_SIM_GT_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_COMPRESS_SIM_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "GTcps.mtx")


def ASSEMBLER_COMPRESS_SIM_B_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_COMPRESS_SIM_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "Bcps.mtx")


def ASSEMBLER_COMPRESS_SIM_C_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_COMPRESS_SIM_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "Ccps.mtx")


def ASSEMBLER_COMPRESS_SIM_INNER_NODES(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_COMPRESS_SIM_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "inner_nodes.dat")


def ASSEMBLER_COMPRESS_SIM_OUTER_NODES(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_COMPRESS_SIM_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "outer_nodes.dat")

ASSEMBLER_COMPRESS_CTRL_DIR = "assembler_compress_ctrl"

def ASSEMBLER_COMPRESS_CTRL_M_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_COMPRESS_CTRL_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "Mcps.mtx")


def ASSEMBLER_COMPRESS_CTRL_MLOWER_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_COMPRESS_CTRL_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "Mlowercps.mtx")


def ASSEMBLER_COMPRESS_CTRL_MUPPER_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_COMPRESS_CTRL_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "Muppercps.mtx")


def ASSEMBLER_COMPRESS_CTRL_S_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_COMPRESS_CTRL_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE),  "Scps.mtx")


def ASSEMBLER_COMPRESS_CTRL_R_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_COMPRESS_CTRL_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "Rcps.mtx")


def ASSEMBLER_COMPRESS_CTRL_K_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_COMPRESS_CTRL_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "Kcps.mtx")


def ASSEMBLER_COMPRESS_CTRL_G_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_COMPRESS_CTRL_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "Gcps.mtx")


def ASSEMBLER_COMPRESS_CTRL_GT_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_COMPRESS_CTRL_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "GTcps.mtx")


def ASSEMBLER_COMPRESS_CTRL_B_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_COMPRESS_CTRL_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "Bcps.mtx")


def ASSEMBLER_COMPRESS_CTRL_C_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_COMPRESS_CTRL_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "Ccps.mtx")


def ASSEMBLER_COMPRESS_CTRL_INNER_NODES(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_COMPRESS_CTRL_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "inner_nodes.dat")


def ASSEMBLER_COMPRESS_CTRL_OUTER_NODES(ref, RE):
    return os.path.join(OUTPUTDIR(), ASSEMBLER_COMPRESS_CTRL_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "outer_nodes.dat")


"""constant for simulation of linearized navier stokes"""
LINEARIZED_SIM_V = "CG"
LINEARIZED_SIM_V_DIM = 2
LINEARIZED_SIM_Q = "CG"
LINEARIZED_SIM_Q_DIM = 1
LINEARIZED_SIM_SAVE_FREQ = 10
LINEARIZED_SIM_INFO = 0.05
LINEARIZED_SIM_DIR = "lqr_sim"

def LINEARIZED_SIM_U_PVD(ref, RE):
    return os.path.join(OUTPUTDIR(), LINEARIZED_SIM_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "u.pvd")

def LINEARIZED_SIM_U_DELTA_PVD(ref, RE):
    return os.path.join(OUTPUTDIR(), LINEARIZED_SIM_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "u_delta.pvd")

def LINEARIZED_SIM_LOG(ref, RE):
    return os.path.join(OUTPUTDIR(), LINEARIZED_SIM_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "log.txt")


"""constants for control of navier stokes"""
LINEARIZED_CTRL_DELTA = - 0.02
LINEARIZED_CTRL_NM_OUTPUT = 1
LINEARIZED_CTRL_NM_RES2_SAVE = 1e-5
LINEARIZED_CTRL_NM_RES2 = 5e-8
LINEARIZED_CTRL_NM_REL2_CHANGE = 1e-10
LINEARIZED_CTRL_NM_REL_CHANGE = 1e-10
LINEARIZED_CTRL_NM_MAXIT = 30
LINEARIZED_CTRL_ADI_OUTPUT = 0
LINEARIZED_CTRL_ADI_RES2 = 1e-15
LINEARIZED_CTRL_ADI_MAXIT = 1000
LINEARIZED_CTRL_V = "CG"
LINEARIZED_CTRL_V_DIM = 2
LINEARIZED_CTRL_Q = "CG"
LINEARIZED_CTRL_Q_DIM = 1
LINEARIZED_CTRL_SAVE_FREQ = 5
LINEARIZED_CTRL_START_CONTROLLING = 0
LINEARIZED_CTRL_INFO = 0.05
LINEARIZED_CTRL_DIR = "lqr_ctrl"

def LINEARIZED_CTRL_U_PVD(ref, RE):
    return os.path.join(OUTPUTDIR(), LINEARIZED_CTRL_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "u.pvd")

def LINEARIZED_CTRL_U_DELTA_PVD(ref, RE):
    return os.path.join(OUTPUTDIR(), LINEARIZED_CTRL_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "u_delta.pvd")

def LINEARIZED_CTRL_Z_CPS_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), LINEARIZED_CTRL_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "Zcps.mtx")

def LINEARIZED_CTRL_KINF_CPS_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), LINEARIZED_CTRL_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "Kinfcps.mtx")

def LINEARIZED_CTRL_LOG(ref, RE):
    return os.path.join(OUTPUTDIR(), LINEARIZED_CTRL_DIR, parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "ctrl.log")


"""constants for bernoulli"""
BERNOULLI_EIGENVALUES = 150
BERNOULLI_MAXIT = 50

def BERNOULLI_FEED0_CPS_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), "bernoulli", parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "Feed0cps.mtx")


"""constant for eigenvalues"""
def EIGEN_SYS_CPS_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), "eigen", parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "eig_sys.mtx")

def EIGEN_RIC_CPS_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), "eigen", parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "eig_ric.mtx")

def EIGEN_BER_CPS_MTX(ref, RE):
    return os.path.join(OUTPUTDIR(), "eigen", parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "eig_ber.mtx")


"""constant for plotter"""
def PLOTTER_LINEARIZED_SIM_LOG1(ref, ending):
    return os.path.join(OUTPUTDIR(), "lqr_sim_plotter", parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "sim1."+ending)

def PLOTTER_LINEARIZED_CTRL_LOG1(ref, ending):
    return os.path.join(OUTPUTDIR(), "lqr_ctrl_plotter", parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "ctrl1."+ending)

def PLOTTER_LINEARIZED_CTRL_LOG2(ref, ending):
    return os.path.join(OUTPUTDIR(), "lqr_ctrl_plotter", parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "ctrl2."+ending)

def PLOTTER_LINEARIZED_CTRL_LOG3(ref, ending):
    return os.path.join(OUTPUTDIR(), "lqr_ctrl_plotter", parameters["refinement_algorithm"],
                        "ref_{0:d}".format(ref), "ctrl3."+ending)



