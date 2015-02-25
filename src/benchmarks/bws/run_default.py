# -*- coding: utf-8 -*-
import src.benchmarks.bws.bws_const as const
from src.benchmarks.aux import *
import sys

if __name__ == "__main__":
    # set Reynoldsnumbers and refinements and Parameters
    REs = range(100, 1000, 50)
    # refs = [1, 2, 3, 4, 5]
    refs = [1]

    const.LINEARIZED_SIM_T = 60.0
    const.LINEARIZED_SIM_T = 1.0
    const.LINEARIZED_SIM_DT = 0.0025
    const.LINEARIZED_CTRL_T = 60.0
    const.LINEARIZED_CTRL_T = 1.0
    const.LINEARIZED_CTRL_DT = 0.0025
    const.LINEARIZED_SIM_INFO = 0.1
    const.LINEARIZED_CTRL_INFO = 0.1
    const.LQR_ADI_OUTPUT = 1

    # read possible inputvar
    if len(sys.argv) == 2:
        lams = [float(sys.argv[1])]
    else:
        lams = [-1.0, -0.5, 0.0, 0.5, 1.0]

    for lam in lams:
        print "lambda ={0:f}".format(lam)
        const.STATIONARY_CONTROL_UPPER.lam = lam
        const.OUTPUTDIR_NAME = "results_bws_lam_{0:3f}".format(lam)

        # build mesh
        print "----------build mesh----------------"
        build_mesh(const, refs)
        print "----------finished build mesh-------"

        # solve stationary
        print "----------newton--------------------"
        solve_newton(const, refs, REs)
        print "----------finished newton-----------"

        # assemble lqr
        REs = [200, 400, 800]
        print "----------assemble------------------"
        assemble_lqr(const, refs, REs)
        print "----------finished assemble---------"

        # simulate
        print "----------simulate------------------"
        simulate(const, refs, REs)
        print "----------finishe simulate----------"

        # solve bernoulli
        print "----------solve bernoulli-----------"
        solve_bernoulli(const, refs, REs)
        print "----------finished solve bernoulli--"

        # solve lqr
        print "----------solve lqr-----------------"
        solve_lqr(const, refs, REs)
        print "----------finishe solve lqr---------"

        # simulate with control
        print "----------start control-------------"
        control(const, refs, REs)
        print "----------finishe control-----------"

        # compute eigenvalues
        # REs = [max(REs)]
        # compute_eigen(const, refs, REs)
	    # break

