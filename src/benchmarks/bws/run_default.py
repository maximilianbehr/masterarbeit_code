# -*- coding: utf-8 -*-
import src.benchmarks.bws.bws_const as const
from src.benchmarks.aux import *
import sys

if __name__ == "__main__":
    # set Reynoldsnumbers and refinements and Parameters
    REs = range(100, 2300, 50)
    refs = [2]
    REscompute = range(min(REs), max(REs), 600)
    const.LINEARIZED_SIM_T = 90.0
    const.LINEARIZED_SIM_DT = 0.002
    const.LINEARIZED_CTRL_T = 90.0
    const.LINEARIZED_CTRL_DT = 0.002
    const.LINEARIZED_SIM_INFO = 0.02
    const.LINEARIZED_CTRL_INFO = 0.02
    const.LQR_ADI_OUTPUT = 1
    const.STATIONARY_LAM = -0.25

    if len(sys.argv) == 4:
        refs = [int(sys.argv[1])]
        const.STATIONARY_LAM = float(sys.argv[2])
        name = sys.argv[3]
        const.OUTPUTDIR_NAME = "results_{0:s}".format(name)
    elif len(sys.argv) == 5:
        refs = [int(sys.argv[1])]
        const.STATIONARY_LAM = float(sys.argv[2])
        name = sys.argv[3]
        const.OUTPUTDIR_NAME = "results_{0:s}".format(name)
        desiredRE = int(sys.argv[4])
        REscompute = [desiredRE]
        REs = range(100, desiredRE, 1000)
        REs.append(desiredRE)


    print "lambda ={0:f}".format(const.STATIONARY_LAM)

    build_mesh(const, refs)
    solve_newton(const, refs, REs)
    REs = REscompute
    assemble_lqr(const, refs, REs)
    compute_condition(const, refs, REs)
    #simulate(const, refs, REs)
    solve_bernoulli(const, refs, REs)
    solve_lqr(const, refs, REs)
    #control(const, refs, REs)
    #compute_eigen(const, refs, REs)