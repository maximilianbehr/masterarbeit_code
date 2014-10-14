from dolfin.cpp.common import info
from dolfin.cpp.mesh import Rectangle,Circle,Mesh,MeshFunction,SubDomain
from dolfin.common.plotting import plot
from dolfin.cpp.io import File
from dolfin.cpp.function import near,between
from dolfin import refine
from math import sqrt

from const.domain import rect, circ
from const.mesh import initial_resolution


# Constants related to the geometry
Rectangledomain     = Rectangle(rect["x0"],rect["y0"],rect["x1"],rect["y1"])
Circledomain        = Circle(circ["x0"],circ["y0"],circ["r"],circ["fragments"])
domain              = Rectangledomain - Circledomain

#Left Boundary
class GammaLeft(SubDomain):
    def inside(self,x,on_boundary):
        return on_boundary and near(x[0],rect["x0"])

#Lower Boundary
class GammaLower(SubDomain):
    def inside(self,x,on_boundary):
        return on_boundary and near(x[1],rect["y0"])

#Right Boundary
class GammaRight(SubDomain):
    def inside(self,x,on_boundary):
        return on_boundary and near(x[0],rect["x1"])

#Upper Boundary
class GammaUpper(SubDomain):
    def inside(self,x,on_boundary):
        return on_boundary and near(x[1],rect["y1"])

#Circle Boundary
class GammaCircle(SubDomain):
    centerx = circ["x0"]
    centery = circ["y0"]
    radius  = circ["r"]
    tresh   = 1.5   #slightly larger use on boundary argument to ensure we get only the boundary of the circle

    def inside(self,x,on_boundary):
        r = sqrt( (x[0]-self.centerx)**2 + (x[1]-self.centery)**2 )
        return on_boundary and r < self.tresh*self.radius

#Circle Controll Boundary 1
class GammaCircleCtrlUpper(SubDomain):
    centerx = circ["x0"]
    centery = circ["y0"]
    radius  = circ["r"]

    uppery  = centery + 6.0/8.0*radius
    lowery  = centery + 1.0/8.0*radius

    def inside(self,x,on_boundary):
        return GammaCircle().inside(x, on_boundary) and between(x[1],(self.lowery,self.uppery)) and self.centerx<x[0]

#Circle Controll Boundary 2
class GammaCircleCtrlLower(SubDomain):
    centerx = circ["x0"]
    centery = circ["y0"]
    radius  = circ["r"]

    uppery  = centery-1.0/8.0*radius
    lowery  = centery-6.0/8.0*radius

    def inside(self,x,on_boundary):
        return GammaCircle().inside(x, on_boundary) and between(x[1],(self.lowery,self.uppery)) and self.centerx<x[0]

#Class for projecting the circle
class CircleProjection(SubDomain):
    centerx = circ["x0"]
    centery = circ["y0"]
    radius  = circ["r"]
    tresh   = 1.5   #slightly larger

    def inside(self,x,on_boundary):
        r       = sqrt( (x[0]-self.centerx)**2 + (x[1]-self.centery)**2 )
        return r < self.tresh*self.radius #dont use on_boundary argument to ensure that a slightly larger circle is transformed

    def snap(self,x):
        r       = sqrt( (x[0]-self.centerx)**2 + (x[1]-self.centery)**2 )
        if r < self.tresh*self.radius:
            x[0]    = self.centerx + (self.radius/r)*(x[0]-self.centerx)
            x[1]    = self.centery + (self.radius/r)*(x[1]-self.centery)


#domain for projection of curved boundary
hole = CircleProjection()


# Mark boundary parts
mesh = Mesh(domain, initial_resolution)
boundaryfunction = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
boundaryfunction.set_all(0)
GammaLeft().mark(boundaryfunction,1)
GammaLower().mark(boundaryfunction,2)
GammaRight().mark(boundaryfunction,3)
GammaUpper().mark(boundaryfunction,4)
GammaCircle().mark(boundaryfunction,5)
GammaCircleCtrlLower().mark(boundaryfunction,6)
GammaCircleCtrlUpper().mark(boundaryfunction,7)

File("mesh.pvd")<<mesh
File("meshboundary.pvd")<<boundaryfunction

##########################################################################

mesh = refine(mesh)
mesh.snap_boundary(hole)
boundaryfunction = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
boundaryfunction.set_all(0)
GammaLeft().mark(boundaryfunction,1)
GammaLower().mark(boundaryfunction,2)
GammaRight().mark(boundaryfunction,3)
GammaUpper().mark(boundaryfunction,4)
GammaCircle().mark(boundaryfunction,5)
GammaCircleCtrlLower().mark(boundaryfunction,6)
GammaCircleCtrlUpper().mark(boundaryfunction,7)

File("ref1mesh.pvd")<<mesh
File("meshboundary1.pvd")<<boundaryfunction

###########################################################################

mesh = refine(mesh)
mesh.snap_boundary(hole)
boundaryfunction = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
boundaryfunction.set_all(0)
GammaLeft().mark(boundaryfunction,1)
GammaLower().mark(boundaryfunction,2)
GammaRight().mark(boundaryfunction,3)
GammaUpper().mark(boundaryfunction,4)
GammaCircle().mark(boundaryfunction,5)
GammaCircleCtrlLower().mark(boundaryfunction,6)
GammaCircleCtrlUpper().mark(boundaryfunction,7)

File("ref2mesh.pvd")<<mesh
File("meshboundary2.pvd")<<boundaryfunction

###########################################################################

mesh = refine(mesh)
mesh.snap_boundary(hole)
boundaryfunction = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
boundaryfunction.set_all(0)
GammaLeft().mark(boundaryfunction,1)
GammaLower().mark(boundaryfunction,2)
GammaRight().mark(boundaryfunction,3)
GammaUpper().mark(boundaryfunction,4)
GammaCircle().mark(boundaryfunction,5)
GammaCircleCtrlLower().mark(boundaryfunction,6)
GammaCircleCtrlUpper().mark(boundaryfunction,7)

File("ref3mesh.pvd")<<mesh
File("meshboundary3.pvd")<<boundaryfunction





















