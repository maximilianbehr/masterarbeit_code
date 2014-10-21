from dolfin import *
from const.nse import *
from const.mesh import *
from mesh import *

import time


def newton_solve(RE,refalg,reflevel):

    # Load mesh and boundary function from file
    mesh                = Mesh(KARMAN_REFN_MESH_FILE(refalg, reflevel))
    boundaryfunction    = MeshFunction("size_t",mesh,KARMAN_REFN_BOUNDARY_FILE(refalg,reflevel))

    # normal vector
    n                   = FacetNormal(mesh)

    # Define function spaces (P2-P1)
    V                   = VectorFunctionSpace(mesh, "CG", 2)
    Q                   = FunctionSpace(mesh, "CG", 1)
    W                   = V*Q

    # define test functions
    (v,q)               = TestFunctions(W)

    # define trial function
    w                   = Function(W)
    (u,p)               = (as_vector((w[0],w[1])),w[2])

    # boundary parts
    ds                  = Measure("ds")[boundaryfunction]

    # Define boundary conditions and rhs
    f                   = Constant(("0","0"))
    u_in                = Expression(("(1-x[1])*x[1]*2","0"))
    noslip_upper        = DirichletBC(W.sub(0), (0, 0), GammaUpper())
    noslip_lower        = DirichletBC(W.sub(0), (0, 0), GammaLower())
    noslip_ball         = DirichletBC(W.sub(0), (0, 0), GammaBall())
    inflow              = DirichletBC(W.sub(0),  u_in , GammaLeft())
    bc                  = [noslip_upper,noslip_lower,noslip_ball,inflow]

    # build weak formulation
    a1                  = inner(grad(u)*u, v)*dx
    a2                  = 1.0/RE*inner(grad(u),grad(v))*dx
    a3                  = -1*p*div(v)*dx
    donothing           = (  -1.0/RE*inner(dot(grad(u),n),v) + inner(p*n,v) )*ds(3)
    cond                = -1*div(u)*q*dx
    rhs                 = inner(f,v)*dx
    F                   = a1 + a2 + a3 + cond

    #
    dw                  = TrialFunction(W)
    dF                  = derivative(F, w, dw)

    nsproblem           = NonlinearVariationalProblem(F, w, bc, dF)
    solver              = NonlinearVariationalSolver(nsproblem)
    solver.parameters["newton_solver"]["maximum_iterations"]=20
    solver.solve()

    (u,p) = w.split()
    File(NSE_VELOCITY_PVD_FILE(RE,refalg,reflevel)) << u
    File(NSE_PRESSURE_PVD_FILE(RE,refalg,reflevel)) << p


if __name__ == "__main__":
    set_log_level(10)
    refalg = "bisection"
    reflevel = 4
    for RE in [1000,2000]:
        print "RE=%d"%RE
        t = time.time()
        newton_solve(RE,refalg,reflevel)
        print "Elapsed =%d s"%(time.time()-t)

