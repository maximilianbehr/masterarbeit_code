# -*- coding: utf-8 -*-
from master.problems.karman.karman import *
from master.problems.karman.karman import __version__
from master.problems.karman.karman import *
from master.problems.karman import karman
from master.problems.karman import karman

DOLFIN_VERSION      = __version__


def STAT_NSE_VELOCITY_PVD_FILE(RE,refalg,reflevel):
    """output for velocity paraview"""
    return "../../results/karman/stat_nse/"+DOLFIN_VERSION+"/%s/ref_%d/RE_%d/velo.pvd"%(refalg,reflevel,RE)

def STAT_NSE_PRESSURE_PVD_FILE(RE,refalg,reflevel):
    """output for pressure paraview"""
    return "../../results/karman/stat_nse/"+DOLFIN_VERSION+"/%s/ref_%d/RE_%d/pressure.pvd"%(refalg,reflevel,RE)


def STAT_NSE_VELOCITY_XML_FILE(RE,refalg,reflevel):
    """output for velocity paraview"""
    return "../../results/karman/stat_nse/"+DOLFIN_VERSION+"/%s/ref_%d/RE_%d/velo.xml"%(refalg,reflevel,RE)

def STAT_NSE_PRESSURE_XML_FILE(RE,refalg,reflevel):
    """output for pressure paraview"""
    return "../../results/karman/stat_nse/"+DOLFIN_VERSION+"/%s/ref_%d/RE_%d/pressure.xml"%(refalg,reflevel,RE)

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
    f                   = Constant((0,0))
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
    cond                = -1*div(u)*q*dx
    rhs                 = inner(f,v)*dx
    F                   = a1 + a2 + a3 + cond + rhs

    #
    dw                  = TrialFunction(W)
    dF                  = derivative(F, w, dw)

    nsproblem           = NonlinearVariationalProblem(F, w, bc, dF)
    solver              = NonlinearVariationalSolver(nsproblem)
    solver.parameters["newton_solver"]["maximum_iterations"]=20
    solver.solve()

    (u,p) = w.split()
    return u,p


if __name__ == "__main__":
    set_log_level(10)
    refalg = "bisection"
    reflevel = 3

    for RE in [100,200]:
        begin("RE=%d"%RE)
        t = time.time()
        u,p=newton_solve(RE,refalg,reflevel)
        print "Elapsed =%d s"%(time.time()-t)
        end()

        try:
                os.remove(STAT_NSE_VELOCITY_PVD_FILE(RE,refalg,reflevel))
                os.remove(STAT_NSE_PRESSURE_PVD_FILE(RE,refalg,reflevel))
                os.remove(STAT_NSE_VELOCITY_XML_FILE(RE,refalg,reflevel))
                os.remove(STAT_NSE_PRESSURE_XML_FILE(RE,refalg,reflevel))
                File(STAT_NSE_VELOCITY_PVD_FILE(RE,refalg,reflevel), "compressed") << u
                File(STAT_NSE_PRESSURE_PVD_FILE(RE,refalg,reflevel), "compressed") << p
                File(STAT_NSE_VELOCITY_XML_FILE(RE,refalg,reflevel)) << u
                File(STAT_NSE_PRESSURE_XML_FILE(RE,refalg,reflevel)) << p
        except OSError:
            pass



