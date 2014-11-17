# -*- coding: utf-8 -*-


# List of solvers
solvers = ["chorin", "css1", "css2", "ipcs", "grpc", "stat_newton"]

# Wrapper for solver classes
def Solver(name, options):
    "Return solver instance for given solver name"
    exec ("from %s import Solver as NamedSolver" % name)
    return NamedSolver(options)
