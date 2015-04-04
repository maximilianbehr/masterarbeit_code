# -*- coding: utf-8 -*-
import src.benchmarks.bws.bws_const as const
from src.benchmarks.aux import *
import sys

if __name__ == "__main__":
    # set Reynoldsnumbers and refinements and Parameters
    REs = range(100, 3000, 50)
    refs = [2]

    const.LINEARIZED_SIM_T = 30.0
    const.LINEARIZED_SIM_DT = 0.002
    const.LINEARIZED_CTRL_T = 30.0
    const.LINEARIZED_CTRL_DT = 0.002
    const.LINEARIZED_SIM_INFO = 0.02
    const.LINEARIZED_CTRL_INFO = 0.02
    const.LQR_ADI_OUTPUT = 1
    lams = [-0.25]

    if len(sys.argv) == 4:
        refs = [int(sys.argv[1])]
        lams = [float(sys.argv[2])]
        name = sys.argv[3]
        const.OUTPUTDIR_NAME = "results_{0:s}".format(name)
    elif len(sys.argv) == 5:
        refs = [int(sys.argv[1])]
        lams = [float(sys.argv[2])]
        name = sys.argv[3]
        const.OUTPUTDIR_NAME = "results_{0:s}".format(name)
        desiredRE = int(sys.argv[4])
        REscompute = [desiredRE]
        REs = range(10, desiredRE, 10)
        REs.append(desiredRE)


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
        #REs = range(min(REs), max(REs), 1000)
        REs = REscompute
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
        #for ref in refs:
        #    for RE in REs:
        #        eig = Eigen(const, ref, RE)
        #        eig.compute_eig_sys()
        #        eig.save()
        #        eig.plot()
