from dolfin.cpp.common import info,begin,end
from dolfin.cpp.mesh import Rectangle,Circle,Mesh,MeshFunction,SubDomain
from dolfin.common.plotting import plot
from dolfin.cpp.io import File
from dolfin.cpp import parameters
from dolfin.cpp.function import near,between
from dolfin import refine
from math import sqrt

from const.domain import rect, circ
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
    boundaryfunction.set_all(0)
    GammaLeft().mark(boundaryfunction,1)
    GammaLower().mark(boundaryfunction,2)
    GammaRight().mark(boundaryfunction,3)
    GammaUpper().mark(boundaryfunction,4)
    GammaBall().mark(boundaryfunction,5)
    GammaBallCtrlLower().mark(boundaryfunction,6)
    GammaBallCtrlUpper().mark(boundaryfunction,7)

    return mesh,boundaryfunction

def refinekarman(mesh):

    #refine given karman
    refmesh = refine(mesh)

    #hole
    hole    = BallProjection()
    refmesh.snap_boundary(hole)

    #mark boundary parts
    boundaryfunction = MeshFunction("size_t", refmesh, refmesh.topology().dim()-1)
    boundaryfunction.set_all(0)
    GammaLeft().mark(boundaryfunction,1)
    GammaLower().mark(boundaryfunction,2)
    GammaRight().mark(boundaryfunction,3)
    GammaUpper().mark(boundaryfunction,4)
    GammaBall().mark(boundaryfunction,5)
    GammaBallCtrlLower().mark(boundaryfunction,6)
    GammaBallCtrlUpper().mark(boundaryfunction,7)

    return refmesh,boundaryfunction



if __name__=="__main__":

    mesh,boundaryfunction = macrokarman()
    File(KARMAN_MACRO_MESH_FILE)<<mesh
    File(KARMAN_MACRO_MESH_PVD_FILE)<<mesh
    File(KARMAN_MACRO_BOUNDARY_FILE)<<boundaryfunction
    File(KARMAN_MACRO_BOUNDARY_PVD_FILE)<<boundaryfunction


    for refalg in parameters.get_range("refinement_algorithm"):
        parameters["refinement_algorithm"]=refalg

        begin("Refinement with %s"%parameters["refinement_algorithm"])
        for ref in range(1,8):
            info("Level %d"%ref)
            mesh, boundaryfunction = refinekarman(mesh)
            File(KARMAN_REFN_MESH_FILE(parameters["refinement_algorithm"],ref))<<mesh
            File(KARMAN_REFN_MESH_PVD_FILE(parameters["refinement_algorithm"],ref))<<mesh
            File(KARMAN_REFN_BOUNDARY_FILE(parameters["refinement_algorithm"],ref))<<boundaryfunction
            File(KARMAN_REFN_BOUNDARY_PVD_FILE(parameters["refinement_algorithm"],ref))<<boundaryfunction

        end("Refinement with %s"%parameters["refinement_algorithm"])





















