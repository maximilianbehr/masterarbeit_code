# -*- coding: utf-8 -*-

import numpy as np
from distutils.version import LooseVersion
from dolfin.cpp.mesh import SubDomain, Mesh, MeshFunction, refine
from dolfin.cpp.function import near, between
from dolfin.cpp.io import File
from dolfin import parameters

import src.karman_const as const


if LooseVersion(const.DOLFIN_VERSION) < LooseVersion("1.4.0"):
    from dolfin.cpp.mesh import Rectangle
    from dolfin.cpp.mesh import Circle
else:
    from mshr import Rectangle
    from mshr import Circle
    from mshr import generate_mesh
    from dolfin import Point


class GammaLeft(SubDomain):
    """Left Boundary"""
    index = const.GAMMA_LEFT_INDICES

    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], const.RECT["x0_0"])


class GammaLower(SubDomain):
    """Lower Boundary"""
    index = const.GAMMA_LOWER_INDICES

    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], const.RECT["x0_1"])


class GammaRight(SubDomain):
    """Right Boundary"""
    index = const.GAMMA_RIGHT_INDICES

    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], const.RECT["x1_0"])


class GammaUpper(SubDomain):
    """Upper Boundary"""
    index = const.GAMMA_UPPER_INDICES

    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], const.RECT["x1_1"])


class GammaBallCtrlUpper(SubDomain):
    """Ball Upper Control Part"""
    index = const.GAMMA_BALL_CTRLUPPER
    upper = const.GAMMA_BALL_CTRLUPPER_UPPER_X1
    lower = const.GAMMA_BALL_CTRLUPPER_LOWER_X1
    thresh = const.GAMMA_BALL_CTRLUPPER_THRESHOLD

    def inside(self, x, on_boundary):
        return (on_boundary and
                GammaBall().inside(x, on_boundary) and
                between(x[1], (self.lower, self.upper))
                and const.CIRCLE["x0"] < x[0])


class GammaBallCtrlLower(SubDomain):
    """Ball Lower Lower Control Part"""
    index = const.GAMMA_BALL_CTRLLOWER
    upper = const.GAMMA_BALL_CTRLLOWER_UPPER_X1
    lower = const.GAMMA_BALL_CTRLLOWER_LOWER_X1
    thresh = const.GAMMA_BALL_CTRLLOWER_THRESHOLD

    def inside(self, x, on_boundary):
        return (on_boundary and
                GammaBall().inside(x, on_boundary) and
                between(x[1], (self.lower, self.upper)) and
                const.CIRCLE["x0"] < x[0])


class GammaBall(SubDomain):
    """Ball"""
    index = const.GAMMA_BALL_INDICES
    thresh = const.GAMMA_BALL_THRESHOLD

    def inside(self, x, on_boundary):
        r = np.sqrt((x[0] - const.CIRCLE["x0"]) ** 2 + (x[1] - const.CIRCLE["x1"]) ** 2)
        return on_boundary and r < self.thresh * const.CIRCLE["r"]


class BallProjection(SubDomain):
    """Ball Projection, special class for projection curved boundary in refinement"""
    thresh = const.GAMMA_BALL_PROJECTION_THRESHOLD

    def inside(self, x, on_boundary):
        r = np.sqrt((x[0] - const.CIRCLE["x0"]) ** 2 + (x[1] - const.CIRCLE["x1"]) ** 2)
        return r < self.thresh * const.CIRCLE["r"]

    def snap(self, x):
        r = np.sqrt((x[0] - const.CIRCLE["x0"]) ** 2 + (x[1] - const.CIRCLE["x1"]) ** 2)
        if r < self.thresh * const.CIRCLE["r"]:
            x[0] = const.CIRCLE["x0"] + (const.CIRCLE["r"] / r) * (x[0] - const.CIRCLE["x0"])
            x[1] = const.CIRCLE["x1"] + (const.CIRCLE["r"] / r) * (x[1] - const.CIRCLE["x1"])


class GammaInner(SubDomain):
    """Inner Edges"""
    index = const.GAMMA_INNER_INDICES


class MeshBuilder():
    """class for refining mesh and meshfunction"""

    def __init__(self):
        print "Dof reordering {0:s}".format(str(parameters["reorder_dofs_serial"]))
        self.mesh = None
        self.boundaryfunction = None
        self.refinelevel = 0

    if LooseVersion(const.DOLFIN_VERSION) < LooseVersion("1.4.0"):
        def _buildmesh(self):
            """use Constants related to the geometry for building CSG"""

            rectangle = Rectangle(const.RECT["x0_0"], const.RECT["x0_1"], const.RECT["x1_0"], const.RECT["x1_1"])
            ball = Circle(const.CIRCLE["x0"], const.CIRCLE["x1"], const.CIRCLE["r"], const.CIRCLE["fragments"])
            domain = rectangle - ball
            self.mesh = Mesh(domain, const.INITIALRESOLUTION)
    else:
        def _buildmesh(self):
            """use Constants related to the geometry for building CSG"""

            p1 = Point(const.RECT["x0_0"], const.RECT["x0_1"])
            p2 = Point(const.RECT["x1_0"], const.RECT["x1_1"])
            p3 = Point(const.CIRCLE["x0"], const.CIRCLE["x1"])

            rectangle = Rectangle(p1, p2)
            ball = Circle(p3, const.CIRCLE["r"])
            domain = rectangle - ball
            self.mesh = generate_mesh(domain, const.INITIALRESOLUTION)


    def _buildboundaryfunction(self):
        """Mark boundary parts"""

        self.boundaryfunction = MeshFunction("size_t", self.mesh, self.mesh.topology().dim() - 1)

        # init all edges with zero
        self.boundaryfunction.set_all(GammaInner().index)

        # mark edges
        GammaLeft().mark(self.boundaryfunction, GammaLeft().index)
        GammaLower().mark(self.boundaryfunction, GammaLower().index)
        GammaRight().mark(self.boundaryfunction, GammaRight().index)
        GammaUpper().mark(self.boundaryfunction, GammaUpper().index)
        GammaBall().mark(self.boundaryfunction, GammaBall().index)
        GammaBallCtrlLower().mark(self.boundaryfunction, GammaBallCtrlLower().index)
        GammaBallCtrlUpper().mark(self.boundaryfunction, GammaBallCtrlUpper().index)


    def refine(self):
        if self.refinelevel == 0:
            self._buildmesh()
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
        File(const.MESH_XML(self.refinelevel), "compressed") << self.mesh
        File(const.MESH_PVD(self.refinelevel), "compressed") << self.mesh
        File(const.MESH_XDMF(self.refinelevel), "compressed") << self.mesh
        File(const.BOUNDARY_XML(self.refinelevel), "compressed") << self.boundaryfunction
        File(const.BOUNDARY_PVD(self.refinelevel), "compressed") << self.boundaryfunction
        File(const.BOUNDARY_XDMF(self.refinelevel), "compressed") << self.boundaryfunction


