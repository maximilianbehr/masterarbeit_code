from __future__ import print_function
from dolfin import *
from karman.problems.const import INSTAT_NSE_CHORIN_PRESSURE_PVD_FILE, INSTAT_NSE_CHORIN_VELOCITY_PVD_FILE
from karman.problems.const import *
from mesh import *

refalg      = "bisection"
reflevel    = 3
RE          = 1


# Print log messages only from the root process in parallel
parameters["std_out_all_processes"] = False;

# Load mesh and boundary function from file
mesh                = Mesh(KARMAN_REFN_MESH_FILE(refalg, reflevel))
boundaryfunction    = MeshFunction("size_t",mesh,KARMAN_REFN_BOUNDARY_FILE(refalg,reflevel))

# Define function spaces (P2-P1)
V = VectorFunctionSpace(mesh, "Lagrange", 2)
Q = FunctionSpace(mesh, "Lagrange", 1)

# Define trial and test functions
u = TrialFunction(V)
p = TrialFunction(Q)
v = TestFunction(V)
q = TestFunction(Q)

# Set parameter values
dt = 0.1
T = 10
nu = 0.01

# Define boundary condition
u_in                = Expression(("(1-x[1])*x[1]*2","0"))
u_noslip            = Constant((0,0))
u_noslip_upper      = DirichletBC(V, u_noslip, GammaUpper())
u_noslip_lower      = DirichletBC(V, u_noslip, GammaLower())
u_noslip_ball       = DirichletBC(V, u_noslip, GammaBall())
u_inflow            = DirichletBC(V, u_in, GammaLeft())
bcu                 = [u_noslip_upper, u_noslip_lower,u_noslip_ball, u_inflow]
#p_in                = Expression("sin(3.0*t)", t=0.0)
#p_inflow            = DirichletBC(Q, p_in,GammaLeft())
#p_outflow           = DirichletBC(Q, 0,  GammaRight())
#bcp                 = [p_inflow,p_outflow]
bcp                 = []


# Create functions
u0 = Function(V)
u1 = Function(V)
p1 = Function(Q)

# Define coefficients
k = Constant(dt)
f = Constant((0, 0))

# Tentative velocity step
F1 = (1/k)*inner(u - u0, v)*dx + inner(grad(u0)*u0, v)*dx + \
     nu*inner(grad(u), grad(v))*dx - inner(f, v)*dx
a1 = lhs(F1)
L1 = rhs(F1)

# Pressure update
a2 = inner(grad(p), grad(q))*dx
L2 = -(1/k)*div(u1)*q*dx

# Velocity update
a3 = inner(u, v)*dx
L3 = inner(u1, v)*dx - k*inner(grad(p1), v)*dx

# Assemble matrices
A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)

# Use amg preconditioner if available
prec = "amg" if has_krylov_solver_preconditioner("amg") else "default"

# files
ufile = File(INSTAT_NSE_CHORIN_VELOCITY_PVD_FILE(1, "bisection", 1))
pfile = File(INSTAT_NSE_CHORIN_PRESSURE_PVD_FILE(1, "bisection", 1))


# Time-stepping
t = dt
while t < T + DOLFIN_EPS:

    # Update pressure boundary condition
    #p_in.t = t

    # Compute tentative velocity step
    begin("Computing tentative velocity")
    b1 = assemble(L1)
    [bc.apply(A1, b1) for bc in bcu]
    solve(A1, u1.vector(), b1, "gmres", "default")
    end()

    # Pressure correction
    begin("Computing pressure correction")
    b2 = assemble(L2)
    if len(bcp) == 0:normalize(b2)
    [bc.apply(A2, b2) for bc in bcp]
    solve(A2, p1.vector(), b2, "cg", prec)
    if len(bcp) == 0:normalize(p1.vector())
    end()

    # Velocity correction
    begin("Computing velocity correction")
    b3 = assemble(L3)
    [bc.apply(A3, b3) for bc in bcu]
    solve(A3, u1.vector(), b3, "gmres", "default")
    end()

    ufile<<(u1,t)
    pfile<<(p1,t)

    # Move to next time step
    u0.assign(u1)
    t += dt
    print("t =", t)

