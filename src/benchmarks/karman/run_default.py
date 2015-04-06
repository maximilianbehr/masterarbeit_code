# -*- coding: utf-8 -*-
from src.benchmarks.aux import *
import src.benchmarks.karman.karman_const as const
import sys

if __name__ == "__main__":
    # set Reynoldsnumbers and refinements and Parameters
    REs = range(10, 200, 10)
    REscompute = range(min(REs), max(REs), 50)
    refs = [2]

    const.LINEARIZED_SIM_T = 20
    const.LINEARIZED_SIM_DT = 0.002
    const.LINEARIZED_CTRL_T = 20
    const.LINEARIZED_CTRL_DT = 0.002
    const.LINEARIZED_SIM_INFO = 0.01
    const.LINEARIZED_CTRL_INFO = 0.1
    const.LQR_ADI_OUTPUT = 1

    if len(sys.argv) == 3:
        refs = [int(sys.argv[1])]
        name = sys.argv[2]
        const.OUTPUTDIR_NAME = "results_{0:s}".format(name)
    elif len(sys.argv) == 4:
        refs = [int(sys.argv[1])]
        name = sys.argv[2]
        const.OUTPUTDIR_NAME = "results_{0:s}".format(name)
        desiredRE = int(sys.argv[3])
        REscompute = [desiredRE]
        REs = range(10, desiredRE, 10)
        REs.append(desiredRE)

    build_mesh(const, refs)
    solve_newton(const, refs, REs)
    REs = REscompute
    assemble_lqr(const, refs, REs)
    simulate(const, refs, REs)
    solve_bernoulli(const, refs, REs)
    solve_lqr(const, refs, REs)
    control(const, refs, REs)

    # compute eigenvalues
    # REs = [max(REs)]
    # compute_eigen(const, refs, REs)



