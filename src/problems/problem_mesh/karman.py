# -*- coding: utf-8 -*-
from numpy import sqrt

from dolfin.cpp.mesh import SubDomain
from dolfin.cpp.mesh import Rectangle
from dolfin.cpp.mesh import Circle
from dolfin.cpp.mesh import refine
from dolfin.cpp.mesh import Mesh
from dolfin.cpp.mesh import MeshFunction
from dolfin.cpp.function import near
from dolfin.cpp.function import between


"constants for rectangular domain"""
rect = {"x0": 0, "x1": 4, "y0": 0, "y1": 1}

"constants for circular domain"
circle = {"x0": 1.5, "y0": 0.5, "r": 0.15, "fragments": 8}

"the resolution of the macro master problem_mesh"
initialresolution = 1


class GammaLeft(SubDomain):
    """Left Boundary"""
    index = 1

    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], rect["x0"])


class GammaLower(SubDomain):
    """Lower Boundary"""
    index = 2

    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], rect["y0"])


class GammaRight(SubDomain):
    """Right Boundary"""
    index = 3

    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], rect["x1"])


class GammaUpper(SubDomain):
    """Upper Boundary"""
    index = 4

    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], rect["y1"])


class GammaBall(SubDomain):
    """Ball"""
    index = 5
    threshold = 1.5  # slightly larger use on boundary argument to ensure we get only the boundary of the Ball

    def inside(self, x, on_boundary):
        r = sqrt((x[0] - circle["x0"]) ** 2 + (x[1] - circle["y0"]) ** 2)
        return on_boundary and self.threshold * circle["r"] > r


class GammaBallCtrlUpper(SubDomain):
    """Ball Upper"""
    index = 6
    uppery = circle["y0"] + 6.0 / 8.0 * circle["r"]
    lowery = circle["y0"] + 1.0 / 8.0 * circle["r"]

    def inside(self, x, on_boundary):
        return GammaBall().inside(x, on_boundary) and between(x[1], (self.lowery, self.uppery)) and circle["x0"] < x[0]


class GammaBallCtrlLower(SubDomain):
    """Ball Lower"""
    index = 7
    uppery = circle["y0"] - 1.0 / 8.0 * circle["r"]
    lowery = circle["y0"] - 6.0 / 8.0 * circle["r"]

    def inside(self, x, on_boundary):
        return GammaBall().inside(x, on_boundary) and between(x[1], (self.lowery, self.uppery)) and circle["x0"] < x[0]


class BallProjection(SubDomain):
    """Ball Projection"""
    threshold = 1.1  # slightly larger

    def inside(self, x, on_boundary):
        r = sqrt((x[0] - circle["x0"]) ** 2 + (x[1] - circle["y0"]) ** 2)
        # dont use on_boundary argument to ensure that a slightly larger Ball is transformed
        return r < self.threshold * circle["r"]

    def snap(self, x):
        r = sqrt((x[0] - circle["x0"]) ** 2 + (x[1] - circle["y0"]) ** 2)
        if r < self.threshold * circle["r"]:
            x[0] = circle["x0"] + (circle["r"] / r) * (x[0] - circle["x0"])
            x[1] = circle["y0"] + (circle["r"] / r) * (x[1] - circle["y0"])


class GammaInner(SubDomain):
    """Inner Edges"""
    index = 0


class MeshBuilder():
    """class for refining mesh and meshfunction"""


    def __init__(self):
        self.mesh = None
        self.boundaryfunction = None
        self.domain = None
        self.refinelevel = 0

        self._buildDomain()
        self._buildMesh()
        self._buildBoundaryFunction()

    def _buildDomain(self):
        """use Constants related to the geometry for building CSG"""
        rectangledomain = Rectangle(rect["x0"], rect["y0"], rect["x1"], rect["y1"])
        balldomain = Circle(circle["x0"], circle["y0"], circle["r"], circle["fragments"])
        self.domain = rectangledomain - balldomain


    def _buildMesh(self):
        self.mesh = Mesh(self.domain, initialresolution)


    def _buildBoundaryFunction(self):
        """Mark boundary parts"""
        self.boundaryfunction = MeshFunction("size_t", self.mesh, self.mesh.topology().dim() - 1)

        # init all edges with zero
        self.boundaryfunction.set_all(GammaInner.index)

        # mark edges
        gammas = [GammaLeft, GammaLower, GammaRight, GammaUpper, GammaBall, GammaBallCtrlLower, GammaBallCtrlUpper]
        for gamma in gammas:
            gamma().mark(self.boundaryfunction, gamma.index)

    def refine(self):
        # refine
        self.mesh = refine(self.mesh)

        # project
        hole = BallProjection()
        self.mesh.snap_boundary(hole)

        self._buildBoundaryFunction()

        #increment refinementlevel
        self.refinelevel += 1




