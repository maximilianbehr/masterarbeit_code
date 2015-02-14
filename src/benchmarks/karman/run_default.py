# -*- coding: utf-8 -*-
import src.benchmarks.karman.karman_const as const
from src.benchmarks.aux import *
from src.aux import print_prof_data

if __name__ == "__main__":
    # set Reynoldsnumbers and refinements and Parameters
    REs = range(10, 120, 10)
    # refs = [1, 2, 3, 4, 5]
    refs = [1]


    #const.LINEARIZED_SIM_DT = 0.005
    #const.LINEARIZED_SIM_T = 10
    const.LINEARIZED_SIM_INFO = 0.1

    # build mesh
    #build_mesh(const, refs)

    # solve stationary
    #solve_newton(const, refs, REs)

    # assemble lqr
    #assemble_lqr(const, refs, REs)

    # simulate
    #simulate(const, refs, REs)

    # solve bernoulli
    #solve_bernoulli(const, refs, REs)

    # solve lqr
    #solve_lqr(const, refs, REs)

    # simulate with control
    #control(const, refs, REs)

    # compute eigenvalues
    compute_eigen(const, refs, REs)

    # plot
    #plot(const, refs, REs)
