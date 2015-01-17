from numpy import linspace
from dolfin import *
from src.mesh.karman import GammaLeft
from src.mesh.karman import GammaLower
from src.mesh.karman import GammaUpper
from src.mesh.karman import GammaBall
from src.mesh.karman import GammaBallCtrlLower
from src.mesh.karman import GammaBallCtrlUpper

import src.karman_const as const


def timestep(T, Umax, nu, mesh, dt=None):
    """Return time step and number of time steps for problem."""
    h = mesh.hmin()

    if dt:
        n = int(T / dt)
        print "Using dt"
    else:
        dt = 0.2 * (h / Umax)
        n = int(T / dt + 1.0)
        dt = T / n
        print "Computing time step according to stability criteria"

    # Compute range
    t_range = linspace(0, T, n + 1)[1:]

    # Report time step
    print 'Number of timesteps:', len(t_range)
    print 'Size of timestep:', dt

    return dt, t_range[0], t_range


def boundary_conditions(V, Q, t):
    # inflow profile
    uin = const.INSTATIONARY_UIN(V, Q, t)

    # noslip at boundary parts
    noslip = Constant((0.0, 0.0))

    # define and collect boundary conditions
    bcu = [DirichletBC(V, uin, GammaLeft()),
           DirichletBC(V, noslip, GammaLower()),
           DirichletBC(V, noslip, GammaUpper()),
           DirichletBC(V, noslip, GammaBall()),
           DirichletBC(V, noslip, GammaBallCtrlLower()),
           DirichletBC(V, noslip, GammaBallCtrlUpper())]

    # define and collect boundary conditions for the pressure
    bcp = []
    return bcu, bcp


def check_divergence(u,V,Q):
    """Check divergence of velocity."""

    # Compute L2 norm of divergence
    print "||div u||_L2 =", norm(u, "Hdiv0")

    # Compute projection of div u into Q_0
    pdivu = project(div(u), Q)
    zero = Constant(0.0)
    bc = DirichletBC(Q, zero, DomainBoundary())
    bc.apply(pdivu.vector())

    # Compute "weak" L2 norm of divergence
    print "||div u||_w  =", sqrt(abs(assemble(pdivu * div(u) * dx)))

