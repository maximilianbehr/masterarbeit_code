from dolfin import *


from src.outputhandler import KarmanOutputHandler
from src.outputhandler import ProblemSolverOutputHandler



kohandler = KarmanOutputHandler()
mesh = Mesh(kohandler.karman_mesh_xml(3))


V = VectorFunctionSpace(mesh,"CG",2)
Q = VectorFunctionSpace(mesh,"CG",1)


dudt = TrialFunction(V)
u = TrialFunction(V)
chi = TrialFunction(Q)



