# -*- coding: utf-8 -*-
from dolfin.cpp.mesh import SubDomain, MeshFunction, refine, cells
from dolfin.cpp.function import near, between
from dolfin.cpp.io import File
from dolfin import parameters
from mshr import Rectangle
from mshr import generate_mesh
from dolfin import Point

import src.benchmarks.bws.bws_const as const


class Gamma1(SubDomain):
    """Left Boundary upper rectangle"""
    index = const.GAMMA1_INDICES

    def inside(self, x, on_boundary):
        return on_boundary and \
               between(x[1], (const.RECTUPPER["x0_1"], const.RECTUPPER["x1_1"])) and \
               near(x[0], const.RECTUPPER["x0_0"])


class Gamma2(SubDomain):
    """lower boundary upper rectangle"""
    index = const.GAMMA2_INDICES

    def inside(self, x, on_boundary):
        return on_boundary and \
               between(x[0], (const.RECTUPPER["x0_0"], const.RECTLOWER["x0_0"])) and \
               near(x[1], const.RECTUPPER["x0_1"])


class Gamma3(SubDomain):
    """upper control boundary lower rectangle left boundary"""
    index = const.GAMMA3_INDICES

    def inside(self, x, on_boundary):
        return on_boundary and \
               near(x[0], const.RECTLOWER["x0_0"]) and \
               between(x[1], (const.RECTUPPER["x0_1"]-const.CONTROLRADIUS, const.RECTUPPER["x0_1"]))


class Gamma4(SubDomain):
    """left boundary lower rectangle"""
    index = const.GAMMA4_INDICES

    def inside(self, x, on_boundary):
        return on_boundary and \
            near(x[0], const.RECTLOWER["x0_0"]) and \
            between(x[1], (const.CONTROLRADIUS, const.RECTUPPER["x0_1"]-const.CONTROLRADIUS))


class Gamma5(SubDomain):
    """lower control boundary lower rectangle left boundary"""
    index = const.GAMMA5_INDICES

    def inside(self, x, on_boundary):
        return on_boundary and \
            near(x[0], const.RECTLOWER["x0_0"]) and \
            between(x[1], (0.0, const.CONTROLRADIUS))


class Gamma6(SubDomain):
    """lower boundary lower rectangle"""
    index = const.GAMMA6_INDICES

    def inside(self, x, on_boundary):
        return on_boundary and \
               between(x[0], (const.RECTLOWER["x0_0"], const.RECTLOWER["x1_0"])) and \
               near(x[1], const.RECTUPPER["x0_0"])


class Gamma7(SubDomain):
    """right boundary both rectangle"""
    index = const.GAMMA7_INDICES

    def inside(self, x, on_boundary):
        return on_boundary and \
        between(x[1], (const.RECTUPPER["x0_0"], const.RECTUPPER["x1_1"])) and \
        near(x[0], const.RECTUPPER["x1_0"])


class Gamma8(SubDomain):
    """upper boundary upper rectangle"""
    index = const.GAMMA8_INDICES

    def inside(self, x, on_boundary):
        return on_boundary and \
        between(x[0], (const.RECTUPPER["x0_0"], const.RECTUPPER["x1_0"])) and \
        near(x[1], const.RECTUPPER["x1_1"])


class GammaShear1(SubDomain):
    """shear layer lower boundary lower rectangle"""
    index = const.GAMMASHEAR1_INDICES

    def inside(self, x, on_boundary):
        return on_boundary and \
               between(x[0], (8*const.MODELHEIGHT, 15*const.MODELHEIGHT)) and \
               near(x[1], 0.0)


class GammaShear2(SubDomain):
    """shear layer lower boundary lower rectangle"""
    index = const.GAMMASHEAR2_INDICES

    def inside(self, x, on_boundary):
        return on_boundary and \
               between(x[0], (7*const.MODELHEIGHT, 15*const.MODELHEIGHT)) and \
               near(x[1], 0.0)


class GammaShear3(SubDomain):
    """shear layer lower boundary lower rectangle"""
    index = const.GAMMASHEAR3_INDICES

    def inside(self, x, on_boundary):
        return on_boundary and \
               between(x[0], (6*const.MODELHEIGHT, 15*const.MODELHEIGHT)) and \
               near(x[1], 0.0)


class GammaInner(SubDomain):
    """Left Boundary upper rectangle"""
    index = const.GAMMA_INNER_INDICES


class MeshBuilder():
    """class for refining mesh and meshfunction"""

    def __init__(self, const):
        print "Dof reordering {0:s}".format(str(parameters["reorder_dofs_serial"]))
        self.mesh = None
        self.boundaryfunction = None
        self.refinelevel = 0
        self.const = const

    def _buildmesh(self):
        """use Constants related to the geometry for building CSG"""

        p1 = Point(self.const.RECTLOWER["x0_0"], self.const.RECTLOWER["x0_1"])
        p2 = Point(self.const.RECTLOWER["x1_0"], self.const.RECTLOWER["x1_1"])
        p3 = Point(self.const.RECTUPPER["x0_0"], self.const.RECTUPPER["x0_1"])
        p4 = Point(self.const.RECTUPPER["x1_0"], self.const.RECTUPPER["x1_1"])
        rectlower = Rectangle(p1, p2)
        rectupper = Rectangle(p3, p4)
        domain = rectlower + rectupper
        self.mesh = generate_mesh(domain, self.const.INITIALRESOLUTION)

    def _buildboundaryfunction(self):
        """Mark boundary parts"""

        self.boundaryfunction = MeshFunction("size_t", self.mesh, self.mesh.topology().dim() - 1)
        self.shear1function = MeshFunction("size_t", self.mesh, self.mesh.topology().dim() - 1)
        self.shear2function = MeshFunction("size_t", self.mesh, self.mesh.topology().dim() - 1)
        self.shear3function = MeshFunction("size_t", self.mesh, self.mesh.topology().dim() - 1)

        # init all edges with zero
        self.boundaryfunction.set_all(GammaInner().index)
        self.shear1function.set_all(GammaInner().index)
        self.shear2function.set_all(GammaInner().index)
        self.shear3function.set_all(GammaInner().index)

        # mark edges for boundary function
        Gamma1().mark(self.boundaryfunction, Gamma1().index)
        Gamma2().mark(self.boundaryfunction, Gamma2().index)
        Gamma4().mark(self.boundaryfunction, Gamma4().index)
        Gamma3().mark(self.boundaryfunction, Gamma3().index)
        Gamma5().mark(self.boundaryfunction, Gamma5().index)
        Gamma6().mark(self.boundaryfunction, Gamma6().index)
        Gamma7().mark(self.boundaryfunction, Gamma7().index)
        Gamma8().mark(self.boundaryfunction, Gamma8().index)

        # mark edges for shear function
        GammaShear1().mark(self.shear1function, GammaShear1().index)
        GammaShear2().mark(self.shear2function, GammaShear2().index)
        GammaShear3().mark(self.shear3function, GammaShear3().index)


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

        # save shear mesh functions
        File(self.const.SHEAR_XML(self.refinelevel, 1), "compressed") << self.shear1function
        File(self.const.SHEAR_PVD(self.refinelevel, 1), "compressed") << self.shear1function
        # File(self.const.SHEAR_XDMF(self.refinelevel, 1), "compressed") << self.shear1function

        File(self.const.SHEAR_XML(self.refinelevel, 2), "compressed") << self.shear2function
        File(self.const.SHEAR_PVD(self.refinelevel, 2), "compressed") << self.shear2function
        # File(self.const.SHEAR_XDMF(self.refinelevel, 2), "compressed") << self.shear2function

        File(self.const.SHEAR_XML(self.refinelevel, 3), "compressed") << self.shear3function
        File(self.const.SHEAR_PVD(self.refinelevel, 3), "compressed") << self.shear3function
        # File(self.const.SHEAR_XDMF(self.refinelevel, 3), "compressed") << self.shear3function


