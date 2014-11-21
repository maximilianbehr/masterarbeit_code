from dolfin import *


from src.outputhandler import KarmanOutputHandler
from src.outputhandler import ProblemSolverOutputHandler
from src.problems.problem_mesh.karman import circle
from src.problems.problem_mesh.karman import GammaBallCtrlLower
from src.problems.problem_mesh.karman import GammaBallCtrlUpper



ref = 4
RE=100
nu = 1.0/RE
penalty_eps = 0.01


kohandler = KarmanOutputHandler()
psohandler = ProblemSolverOutputHandler("karman", "stat_newton")

mesh = Mesh(kohandler.karman_mesh_xml(ref))
boundaryfunction = MeshFunction("size_t",mesh,kohandler.karman_boundary_xml(ref))



g = Expression(("1/r * (x[0]-x0)", "1/r * (x[1]-y0)"), r = circle["r"], x0 = circle["x0"], y0 = circle["y0"] )

V = VectorFunctionSpace(mesh,"CG",2)
Q = FunctionSpace(mesh,"CG",1)

u_stat = Function(V,psohandler.u_xml(ref,RE))


#trial functions
dudt = TrialFunction(V)
u = TrialFunction(V)
p = TrialFunction(Q)


#test functions
w_test = TestFunction(V)
p_test = TestFunction(Q)


#weak formulation of linearized navier stokes eqn
ds = Measure('ds')[boundaryfunction]

M_var = inner(dudt, w_test)*dx
S_var = inner(grad(u),grad(w_test))*dx #1/nu spaeter ranmultiplizieren

B_upper_var = inner(g,w_test)*ds(GammaBallCtrlUpper.index)#1/nu *1/eps spaeter ran
B_lower_var = inner(g,w_test)*ds(GammaBallCtrlLower.index)

M_upper_var = inner(u,w_test)*ds(GammaBallCtrlUpper.index)#1/nu *1/eps spaeter ran
M_lower_var = inner(u,w_test)*ds(GammaBallCtrlLower.index)

K_var = inner(grad(u_stat)*u,w_test)*dx

R_var = inner(grad(u)*u_stat,w_test)*dx

G_var = -1*p*div(w_test)*dx


Gt_var = div(u)*p_test*dx


#assemble matrices

M = assemble(M_var)
S = assemble(S_var)
B_upper = assemble(B_upper_var)
B_lower = assemble(B_lower_var)
M_upper = assemble(M_upper_var)
M_lower = assemble(M_lower_var)
K = assemble(K_var)
R = assemble(R_var)
G = assemble(G_var)
Gt = assemble(Gt_var)