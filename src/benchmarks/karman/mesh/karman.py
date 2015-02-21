# -*- coding: utf-8 -*-

import numpy as np
from dolfin.cpp.mesh import SubDomain, Mesh, MeshFunction, refine, cells
from dolfin.cpp.function import near, between
from dolfin.cpp.io import File
from dolfin import parameters

from src.benchmarks.karman.karman_const import RECT
from src.benchmarks.karman.karman_const import CIRCLE
from src.benchmarks.karman.karman_const import GAMMA_BALL_CTRLUPPER_UPPER_X1
from src.benchmarks.karman.karman_const import GAMMA_BALL_CTRLUPPER_LOWER_X1
from src.benchmarks.karman.karman_const import GAMMA_BALL_CTRLUPPER_THRESHOLD
from src.benchmarks.karman.karman_const import GAMMA_BALL_CTRLLOWER_UPPER_X1
from src.benchmarks.karman.karman_const import GAMMA_BALL_CTRLLOWER_LOWER_X1
from src.benchmarks.karman.karman_const import GAMMA_BALL_CTRLLOWER_THRESHOLD
from src.benchmarks.karman.karman_const import GAMMA_BALL_THRESHOLD
from src.benchmarks.karman.karman_const import GAMMA_BALL_PROJECTION_THRESHOLD
from mshr import Rectangle
from mshr import Circle
from mshr import generate_mesh
from dolfin import Point


class GammaLeft(SubDomain):
    """Left Boundary"""

    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], RECT["x0_0"])


class GammaLower(SubDomain):
    """Lower Boundary"""

    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], RECT["x0_1"])


class GammaRight(SubDomain):
    """Right Boundary"""

    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], RECT["x1_0"])


class GammaUpper(SubDomain):
    """Upper Boundary"""

    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], RECT["x1_1"])


class GammaBallCtrlUpper(SubDomain):
    """Ball Upper Control Part"""
    upper = GAMMA_BALL_CTRLUPPER_UPPER_X1
    lower = GAMMA_BALL_CTRLUPPER_LOWER_X1
    thresh = GAMMA_BALL_CTRLUPPER_THRESHOLD

    def inside(self, x, on_boundary):
        return (on_boundary and
                GammaBall().inside(x, on_boundary) and
                between(x[1], (self.lower, self.upper))
                and CIRCLE["x0"] < x[0])


class GammaBallCtrlLower(SubDomain):
    """Ball Lower Lower Control Part"""
    upper = GAMMA_BALL_CTRLLOWER_UPPER_X1
    lower = GAMMA_BALL_CTRLLOWER_LOWER_X1
    thresh = GAMMA_BALL_CTRLLOWER_THRESHOLD

    def inside(self, x, on_boundary):
        return (on_boundary and
                GammaBall().inside(x, on_boundary) and
                between(x[1], (self.lower, self.upper)) and
                CIRCLE["x0"] < x[0])


class GammaBall(SubDomain):
    """Ball"""
    thresh = GAMMA_BALL_THRESHOLD

    def inside(self, x, on_boundary):
        r = np.sqrt((x[0] - CIRCLE["x0"]) ** 2 + (x[1] - CIRCLE["x1"]) ** 2)
        return on_boundary and r < self.thresh * CIRCLE["r"]


class BallProjection(SubDomain):
    """Ball Projection, special class for projection curved boundary in refinement"""
    thresh = GAMMA_BALL_PROJECTION_THRESHOLD

    def inside(self, x, on_boundary):
        r = np.sqrt((x[0] - CIRCLE["x0"]) ** 2 + (x[1] - CIRCLE["x1"]) ** 2)
        return r < self.thresh * CIRCLE["r"]

    def snap(self, x):
        r = np.sqrt((x[0] - CIRCLE["x0"]) ** 2 + (x[1] - CIRCLE["x1"]) ** 2)
        if r < self.thresh * CIRCLE["r"]:
            x[0] = CIRCLE["x0"] + (CIRCLE["r"] / r) * (x[0] - CIRCLE["x0"])
            x[1] = CIRCLE["x1"] + (CIRCLE["r"] / r) * (x[1] - CIRCLE["x1"])



class MeshBuilder():
    """class for refining mesh and meshfunction"""

    def __init__(self, const):
        print "Dof reordering {0:s}".format(str(parameters["reorder_dofs_serial"]))
        self.mesh =    None
        self.boundaryfunction = None
        self.refinelevel = 0
        self.const = const

    def _buildmesh(self):
        """use Constants related to the geometry for building CSG"""
        p1 = Point(self.const.RECT["x0_0"], self.const.RECT["x0_1"])
        p2 = Point(self.const.RECT["x1_0"], self.const.RECT["x1_1"])
        p3 = Point(self.const.CIRCLE["x0"], self.const.CIRCLE["x1"])

        rectangle = Rectangle(p1, p2)
        ball = Circle(p3, self.const.CIRCLE["r"])
        domain = rectangle - ball
        self.mesh = generate_mesh(domain, self.const.INITIALRESOLUTION)


    def _buildboundaryfunction(self):
        """Mark boundary parts"""

        self.boundaryfunction = MeshFunction("size_t", self.mesh, self.mesh.topology().dim() - 1)

        # init all edges with zero
        self.boundaryfunction.set_all(self.const.GAMMA_INNER_INDICES)

        # mark edges
        GammaLeft().mark(self.boundaryfunction, self.const.GAMMA_LEFT_INDICES)
        GammaLower().mark(self.boundaryfunction, self.const.GAMMA_LOWER_INDICES)
        GammaRight().mark(self.boundaryfunction, self.const.GAMMA_RIGHT_INDICES)
        GammaUpper().mark(self.boundaryfunction, self.const.GAMMA_UPPER_INDICES)
        GammaBall().mark(self.boundaryfunction, self.const.GAMMA_BALL_INDICES)
        GammaBallCtrlLower().mark(self.boundaryfunction, self.const.GAMMA_BALL_CTRLLOWER_INDICES)
        GammaBallCtrlUpper().mark(self.boundaryfunction, self.const.GAMMA_BALL_CTRLUPPER_INDICES)


    def refine(self):
        if self.refinelevel == 0:
            self._buildmesh()

            # local refinement
            for i in range(self.const.LOCALREFINEMENTS):
                cellmarkers = MeshFunction("bool", self.mesh, self.mesh.topology().dim())
                cellmarkers.set_all(False)
                for c in cells(self.mesh):
                    p = c.midpoint()
                    if self.const.LOCALREFINE(p):
                        cellmarkers[c] = True
                self.mesh = refine(self.mesh, cellmarkers)

            self._buildboundaryfunction()
        else:
            # refine
            self.mesh = refine(self.mesh)

            # project
            hole = BallProjection()
            self.mesh.snap_boundary(hole)

            # build boundary function
            self._buildboundaryfunction()

        # increment refinementlevel
        self.refinelevel += 1

    def save(self):
        File(self.const.MESH_XML(self.refinelevel), "compressed") << self.mesh
        File(self.const.MESH_PVD(self.refinelevel), "compressed") << self.mesh
        # File(self.const.MESH_XDMF(self.refinelevel), "compressed") << self.mesh
        File(self.const.BOUNDARY_XML(self.refinelevel), "compressed") << self.boundaryfunction
        File(self.const.BOUNDARY_PVD(self.refinelevel), "compressed") << self.boundaryfunction
        # File(self.const.BOUNDARY_XDMF(self.refinelevel), "compressed") << self.boundaryfunction


