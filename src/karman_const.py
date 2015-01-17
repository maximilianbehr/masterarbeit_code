# -*- coding: utf-8 -*-
import dolfin
from dolfin import parameters
from dolfin import Expression
from dolfin import Constant
import os
import socket


"""dolfin version"""
DOLFIN_VERSION = dolfin.__version__


"""turn dof reordering off"""
parameters["reorder_dofs_serial"] = False

"""set refinement algorithm"""
parameters["refinement_algorithm"] = "recursive_bisection"

"""Optimzation parameters"""
parameters["form_compiler"]["optimize"] = True
parameters["form_compiler"]["cpp_optimize"] = True
parameters["krylov_solver"]["absolute_tolerance"] = 1e-25
parameters["krylov_solver"]["relative_tolerance"] = 1e-14
parameters["krylov_solver"]["monitor_convergence"] = False


"""rectangular domain"""
RECT = {"x0_0": 0, "x1_0": 5.0, "x0_1": 0, "x1_1": 1.0}

"""circular obstacle domain"""
CIRCLE = {"x0": 0.5, "x1": 0.5, "r": 0.15, "fragments": 8}

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
GAMMA_BALL_CTRLLOWER = 6
GAMMA_BALL_CTRLUPPER = 7


"""Threshold for projection of the boundary of the obstacle"""
GAMMA_BALL_PROJECTION_THRESHOLD = 1.1

assert GAMMA_BALL_PROJECTION_THRESHOLD > 1, "Constant must be slightly larger than 1"


"""Output directory"""
def OUTPUTDIR():
    dirname = os.path.join("results", DOLFIN_VERSION)
    host = socket.gethostname()
    if host in ["pc747", "pc633", "pc800"]:
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
STATIONARY_P0 = Constant(0.0)
STATIONARY_UIN = Expression(("4*(1-x[1])*x[1]", "0.0"))
STATIONARY_NEWTON_STEPS = 40
STATIONARY_NEWTON_ABS_TOL = 1e-12
STATIONARY_NEWTON_REL_TOL = 1e-14


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

























