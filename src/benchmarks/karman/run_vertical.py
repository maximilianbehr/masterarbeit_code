# -*- coding: utf-8 -*-
import src.benchmarks.karman.karman_const as const
from src.benchmarks.aux import *
from src.aux import print_prof_data
import sys

if __name__ == "__main__":
    # set Reynoldsnumbers and refinemnts
    REs = range(10, 500, 10)
    refs = [1]

    const.LINEARIZED_SIM_T = 20
    const.LINEARIZED_SIM_DT = 0.002
    const.LINEARIZED_CTRL_T = 20
    const.LINEARIZED_CTRL_DT = 0.002
    const.LINEARIZED_SIM_INFO = 0.02
    const.LINEARIZED_CTRL_INFO = 0.02
    const.LQR_ADI_OUTPUT = 1
    const.ASSEMBLER_OBSERVER_POINTS = [(3.5, 0.5), (4.0, 0.5)]

    if len(sys.argv) == 3:
        refs = [int(sys.argv[1])]
        name = sys.argv[2]
        const.OUTPUTDIR_NAME = "results_{0:s}".format(name)

    build_mesh(const, refs)
    solve_newton(const, refs, REs)
    REs = range(min(REs), max(REs), 70)
    assemble_lqr(const, refs, REs)
    simulate(const, refs, REs)
    solve_bernoulli(const, refs, REs)
    solve_lqr(const, refs, REs)
    control(const, refs, REs)
    # compute_eigen(const, refs, [max(REs)])
    # plot(const, refs, REs)
