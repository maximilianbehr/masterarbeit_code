# -*- coding: utf-8 -*-

from distutils.version import LooseVersion
from dolfin.cpp.mesh import SubDomain, Mesh, MeshFunction, refine
from dolfin.cpp.function import near, between
from dolfin.cpp.io import File
from dolfin import parameters

import src.bws.bws_const as const


if LooseVersion(const.DOLFIN_VERSION) < LooseVersion("1.4.0"):
    from dolfin.cpp.mesh import Rectangle
else:
    from mshr import Rectangle
    from mshr import generate_mesh
    from dolfin import Point


class Gamma1(SubDomain):
    """Left Boundary upper rectangle"""
    index = const.GAMMA1_INDICES

    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], const.RECTUPPER["x0_0"])


class Gamma2(SubDomain):
    """lower boundary upper rectangle"""
    index = const.GAMMA2_INDICES

    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], const.RECTUPPER["x0_1"])


class Gamma3(SubDomain):
    """left boundary lower rectangle"""
    index = const.GAMMA3_INDICES

    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], const.RECTLOWER["x0_0"])


class Gamma4(SubDomain):
    """lower boundary lower rectangle"""
    index = const.GAMMA4_INDICES

    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], const.RECTLOWER["x0_1"])


class Gamma5(SubDomain):
    """right boundary lower and upper rectangle"""
    index = const.GAMMA5_INDICES

    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], const.RECTLOWER["x1_0"])


class Gamma6(SubDomain):
    """upper boundary upper rectangle"""
    index = const.GAMMA6_INDICES

    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], const.RECTUPPER["x1_1"])


class GammaInner(SubDomain):
    """Left Boundary upper rectangle"""
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
            rectlower = Rectangle(const.RECTLOWER["x0_0"], const.RECTLOWER["x0_1"], const.RECTLOWER["x1_0"], const.RECTLOWER["x1_1"])
            rectupper = Rectangle(const.RECTUPPER["x0_0"], const.RECTUPPER["x0_1"], const.RECTUPPER["x1_0"], const.RECTUPPER["x1_1"])
            domain = rectlower + rectupper
            self.mesh = Mesh(domain, const.INITIALRESOLUTION)
    else:
        def _buildmesh(self):
            """use Constants related to the geometry for building CSG"""

            p1 = Point(const.RECTLOWER["x0_0"], const.RECTLOWER["x0_1"])
            p2 = Point(const.RECTLOWER["x1_0"], const.RECTLOWER["x1_1"])
            p3 = Point(const.RECTUPPER["x0_0"], const.RECTUPPER["x0_1"])
            p4 = Point(const.RECTUPPER["x1_0"], const.RECTUPPER["x1_1"])
            rectlower = Rectangle(p1, p2)
            rectupper = Rectangle(p3, p4)
            domain = rectlower + rectupper
            self.mesh = generate_mesh(domain, const.INITIALRESOLUTION)

    def _buildboundaryfunction(self):
        """Mark boundary parts"""

        self.boundaryfunction = MeshFunction("size_t", self.mesh, self.mesh.topology().dim() - 1)

        # init all edges with zero
        self.boundaryfunction.set_all(GammaInner().index)

        # mark edges
        Gamma1().mark(self.boundaryfunction, Gamma1().index)
        Gamma2().mark(self.boundaryfunction, Gamma2().index)
        Gamma3().mark(self.boundaryfunction, Gamma3().index)
        Gamma4().mark(self.boundaryfunction, Gamma4().index)
        Gamma5().mark(self.boundaryfunction, Gamma5().index)
        Gamma6().mark(self.boundaryfunction, Gamma6().index)

    def refine(self):
        if self.refinelevel == 0:
            self._buildmesh()
            self._buildboundaryfunction()
        else:
            # refine
            self.mesh = refine(self.mesh)

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


