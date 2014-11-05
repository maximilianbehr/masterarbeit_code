from dolfin import *
from karman.problems.const import INSTAT_NSE_CSS_PRESSURE_PVD_FILE, INSTAT_NSE_CSS_VELOCITY_PVD_FILE, epsilon, sigma
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


#Define function spaces
V = VectorFunctionSpace(mesh, "CG", 2)
Q = FunctionSpace(mesh, "CG", 1)


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
bcpsi               = homogenize(bcp)


# Test and trial functions
v = TestFunction(V)
q = TestFunction(Q)
u = TrialFunction(V)
p = TrialFunction(Q)

# Set parameter values
dt  = 0.01
T   = 10
nu  = 0.1


# Functions

u0  = interpolate(Constant((0,0)),V)
u1  = interpolate(Constant((0,0)),V)
p0  = interpolate(Constant(0),Q)
p1  = interpolate(Constant(0),Q)
p2  = interpolate(Constant(0),Q)
k   = Constant(dt)
n   = FacetNormal(mesh)
f   = Constant((0,0))
psi = Function(Q)
pbar = Constant(0)

#u0  = interpolate(u0, V)
#u1  = interpolate(u0, V)
#p0  = interpolate(p0, Q)
#p1  = interpolate(p0, Q)
#p2  = interpolate(p0, Q)
#nu  = Constant(problem.nu)
#k   = Constant(dt)
#n   = FacetNormal(mesh)
#f   = Constant(0,0)
#psi = Function(Q)

# Tentative pressure
#if self.order == 1:
#    ps = p1
#else:
#    ps = 2*p1 - p0
#take first order
ps = p1
#ps = 2*p1 - p0


# Tentative velocity step
F1 = (1/k)*inner(v, u - u0)*dx + inner(v, grad(u0)*u0)*dx \
    + inner(epsilon(v), sigma(u, ps, nu))*dx \
    + inner(v, pbar*n)*ds \
    - inner(v, f)*dx
a1 = lhs(F1)
L1 = rhs(F1)

# Pressure correction
a2 = inner(grad(q), grad(p))*dx
L2 = (1/k)*inner(grad(q), u1 - u0)*dx - (1/k)*inner(q*n, u1 - u0)*ds

# Pressure update
a3 = q*p*dx
L3 = q*(ps + psi - nu*div(u1))*dx

# Assemble matrices
A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)


# Use amg preconditioner if available
prec = "amg" if has_krylov_solver_preconditioner("amg") else "default"


# Time loop
t=0;
while t < T + DOLFIN_EPS:

    # Compute tentative velocity step
    b = assemble(L1)
    [bc.apply(A1, b) for bc in bcu]
    solve(A1, u1.vector(), b, "gmres", "ilu")

    # Compute pressure correction
    b = assemble(L2)
    if len(bcp) == 0: normalize(b)
    [bc.apply(A2, b) for bc in bcpsi]
    solve(A2, psi.vector(), b, "gmres", prec)
    if len(bcp) == 0 : normalize(psi.vector())

    # Compute updated pressure
    b = assemble(L3)
    if len(bcp) == 0: normalize(b)
    [bc.apply(A3, b) for bc in bcp]
    solve(A3, p2.vector(), b, "gmres", "ilu")

    # Update
    u0.assign(u1)
    p0.assign(p1)
    p1.assign(p2)

    print "t =%f" %t
    t+=dt


#return u1, p2






