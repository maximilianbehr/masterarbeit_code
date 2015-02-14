# -*- coding: utf-8 -*-
import src.benchmarks.karman.karman_const as const
from src.benchmarks.aux import *
from src.aux import print_prof_data


if __name__ == "__main__":

    #refs = [2]
    #REs = range(40, 200, 20)
    refs = [1]
    REs = range(60, 110, 20)


    const.OUTPUTDIR_NAME = "results_bernoulli"
    const.LINEARIZED_SIM_T = 20
    const.LINEARIZED_SIM_DT = 0.005
    const.LINEARIZED_CTRL_T = 20
    const.LINEARIZED_CTRL_DT = 0.001
    const.LINEARIZED_SIM_PERTUBATIONEPS = 0.25
    const.LINEARIZED_CTRL_PERTUBATIONEPS = 0.25
    # const.ASSEMBLER_OBSERVER_POINTS = [(2.5, 0.5), (3.0, 0.5)]
    const.BERNOULLI_FEED0_CPS_MTX = lambda ref,RE: const.ASSEMBLER_COMPRESS_CTRL_NAME_MTX(ref, "Kinf", RE)

    build_mesh(const, refs)
    solve_newton(const, refs, REs)
    assemble_lqr(const, refs, REs)
    simulate(const, refs, REs)
    solve_bernoulli(const, refs, REs)
    control(const, refs, REs)
    # compute_eigen(const, refs, REs)
    # plot(const, refs, REs)
