# -*- coding: utf-8 -*-
from dolfin.cpp.common import info,begin,end
from dolfin.cpp.mesh import Rectangle,Circle,Mesh,MeshFunction,SubDomain
from dolfin.cpp.common import *
from dolfin.common.plotting import plot
from dolfin.cpp.io import File
from dolfin.cpp import parameters
from dolfin.cpp.function import near,between
from dolfin import refine
from math import sqrt

from const.mesh import *


class GammaLeft(SubDomain):
    """Left Boundary"""
    def inside(self,x,on_boundary):
        return on_boundary and near(x[0],rect["x0"])

class GammaLower(SubDomain):
    """Lower Boundary"""
    def inside(self,x,on_boundary):
        return on_boundary and near(x[1],rect["y0"])

class GammaRight(SubDomain):
    """Right Boundary"""
    def inside(self,x,on_boundary):
        return on_boundary and near(x[0],rect["x1"])

class GammaUpper(SubDomain):
    """Upper Boundary"""
    def inside(self,x,on_boundary):
        return on_boundary and near(x[1],rect["y1"])

class GammaBall(SubDomain):
    """Ball"""
    centerx = circ["x0"]
    centery = circ["y0"]
    radius  = circ["r"]
    tresh   = 1.5   #slightly larger use on boundary argument to ensure we get only the boundary of the Ball

    def inside(self,x,on_boundary):
        r = sqrt( (x[0]-self.centerx)**2 + (x[1]-self.centery)**2 )
        return on_boundary and r < self.tresh*self.radius

class GammaBallCtrlUpper(SubDomain):
    """Ball Upper"""
    centerx = circ["x0"]
    centery = circ["y0"]
    radius  = circ["r"]

    uppery  = centery + 6.0/8.0*radius
    lowery  = centery + 1.0/8.0*radius

    def inside(self,x,on_boundary):
        return GammaBall().inside(x, on_boundary) and between(x[1],(self.lowery,self.uppery)) and self.centerx<x[0]

class GammaBallCtrlLower(SubDomain):
    """Ball Lower"""
    centerx = circ["x0"]
    centery = circ["y0"]
    radius  = circ["r"]

    uppery  = centery-1.0/8.0*radius
    lowery  = centery-6.0/8.0*radius

    def inside(self,x,on_boundary):
        return GammaBall().inside(x, on_boundary) and between(x[1],(self.lowery,self.uppery)) and self.centerx<x[0]

class BallProjection(SubDomain):
    """Ball Projection"""
    centerx = circ["x0"]
    centery = circ["y0"]
    radius  = circ["r"]
    tresh   = 1.5   #slightly larger

    def inside(self,x,on_boundary):
        r       = sqrt( (x[0]-self.centerx)**2 + (x[1]-self.centery)**2 )
        return r < self.tresh*self.radius #dont use on_boundary argument to ensure that a slightly larger Ball is transformed

    def snap(self,x):
        r       = sqrt( (x[0]-self.centerx)**2 + (x[1]-self.centery)**2 )
        if r < self.tresh*self.radius:
            x[0]    = self.centerx + (self.radius/r)*(x[0]-self.centerx)
            x[1]    = self.centery + (self.radius/r)*(x[1]-self.centery)

def macrokarman():

    # use Constants related to the geometry for building CSG
    Rectangledomain     = Rectangle(rect["x0"],rect["y0"],rect["x1"],rect["y1"])
    Balldomain          = Circle(circ["x0"],circ["y0"],circ["r"],circ["fragments"])
    domain              = Rectangledomain - Balldomain

    # Mark boundary parts
    mesh                = Mesh(domain, initial_resolution)
    boundaryfunction    = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
    boundaryfunction.set_all(GAMMAINNER)
    GammaLeft().mark(boundaryfunction,GAMMALEFT)
    GammaLower().mark(boundaryfunction,GAMMALOWER)
    GammaRight().mark(boundaryfunction,GAMMARIGHT)
    GammaUpper().mark(boundaryfunction,GAMMAUPPER)
    GammaBall().mark(boundaryfunction,GAMMABALL)
    GammaBallCtrlLower().mark(boundaryfunction,GAMMABALLCTRLUPPER)
    GammaBallCtrlUpper().mark(boundaryfunction,GAMMABALLCTRLLOWER)

    return mesh,boundaryfunction

def refinekarman(mesh):

    #refine given karman
    refmesh = refine(mesh)

    #hole
    hole    = BallProjection()
    refmesh.snap_boundary(hole)

    #mark boundary parts
    boundaryfunction = MeshFunction("size_t", refmesh, refmesh.topology().dim()-1)
    boundaryfunction.set_all(GAMMAINNER)
    GammaLeft().mark(boundaryfunction,GAMMALEFT)
    GammaLower().mark(boundaryfunction,GAMMALOWER)
    GammaRight().mark(boundaryfunction,GAMMARIGHT)
    GammaUpper().mark(boundaryfunction,GAMMAUPPER)
    GammaBall().mark(boundaryfunction,GAMMABALL)
    GammaBallCtrlLower().mark(boundaryfunction,GAMMABALLCTRLUPPER)
    GammaBallCtrlUpper().mark(boundaryfunction,GAMMABALLCTRLLOWER)

    return refmesh,boundaryfunction


if __name__=="__main__":

    #for refalg in parameters.get_range("refinement_algorithm"):
    #for refalg in ["bisection", "iterative_bisection", "recursive_bisection", "regular_cut"]:
    for refalg in ["bisection"]:

        parameters["refinement_algorithm"]=refalg

        mesh,boundaryfunction = macrokarman()

        File(KARMAN_MACRO_MESH_FILE,"compressed")<<mesh
        File(KARMAN_MACRO_BOUNDARY_FILE,"compressed")<<boundaryfunction

        File(KARMAN_MACRO_MESH_FILE_XDMF,"compressed")<<mesh
        File(KARMAN_MACRO_BOUNDARY_FILE_XDMF,"compressed")<<boundaryfunction

        File(KARMAN_MACRO_MESH_PVD_FILE,"compressed")<<mesh
        File(KARMAN_MACRO_BOUNDARY_PVD_FILE,"compressed")<<boundaryfunction

        begin("Refinement with %s"%parameters["refinement_algorithm"])
        for ref in range(1,8):
            info("Level %d"%ref)
            mesh, boundaryfunction = refinekarman(mesh)

            File(KARMAN_REFN_MESH_FILE(parameters["refinement_algorithm"],ref),"compressed")<<mesh
            File(KARMAN_REFN_BOUNDARY_FILE(parameters["refinement_algorithm"],ref),"compressed")<<boundaryfunction

            File(KARMAN_REFN_MESH_FILE_XDMF(parameters["refinement_algorithm"],ref),"compressed")<<mesh
            File(KARMAN_REFN_BOUNDARY_FILE_XDMF(parameters["refinement_algorithm"],ref),"compressed")<<boundaryfunction

            File(KARMAN_REFN_MESH_PVD_FILE(parameters["refinement_algorithm"],ref),"compressed")<<mesh
            File(KARMAN_REFN_BOUNDARY_PVD_FILE(parameters["refinement_algorithm"],ref),"compressed")<<boundaryfunction

        end()





















