# -*- coding: utf-8 -*-
import src.benchmarks.bws.bws_const as const
from src.benchmarks.aux import *
import sys

if __name__ == "__main__":
    # set Reynoldsnumbers and refinements and Parameters
    REs = range(100, 5000, 100)
    refs = [1]

    const.LINEARIZED_SIM_T = 40.0
    const.LINEARIZED_SIM_DT = 0.00025
    const.LINEARIZED_CTRL_T = 40.0
    const.LINEARIZED_CTRL_DT = 0.00025
    const.LINEARIZED_SIM_INFO = 0.1
    const.LINEARIZED_CTRL_INFO = 0.1
    const.LQR_ADI_OUTPUT = 1
    lams = [-0.25]

    if len(sys.argv) == 3:
        refs = [int(sys.argv[1])]
        lams = [float(sys.argv[2])]
        name = sys.argv[3]
        const.OUTPUTDIR_NAME = "results_{0:s}".format(name)


    for lam in lams:
        print "lambda ={0:f}".format(lam)
        const.STATIONARY_LAM = lam

        # build mesh
        print "----------build mesh----------------"
        build_mesh(const, refs)
        print "----------finished build mesh-------"

        # solve stationary
        print "----------newton--------------------"
        solve_newton(const, refs, REs)
        print "----------finished newton-----------"

        # assemble lqr
        REs = range(min(REs), max(REs), 600)
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
