# -*- coding: utf-8 -*-
import src.benchmarks.karman.karman_const as const
from src.benchmarks.aux import *
import sys

if __name__ == "__main__":
    # set Reynoldsnumbers and refinemnts
    REs = range(10, 200, 10)
    refs = [1]

    const.LINEARIZED_SIM_T = 20
    const.LINEARIZED_SIM_DT = 0.002
    const.LINEARIZED_CTRL_T = 20
    const.LINEARIZED_CTRL_DT = 0.002
    const.LINEARIZED_SIM_INFO = 0.02
    const.LINEARIZED_CTRL_INFO = 0.02
    const.LQR_ADI_OUTPUT = 1
    # const.ASSEMBLER_OBSERVER_POINTS = [(3.5, 0.5), (4.0, 0.5)]
    const.OUTPUTDIR_NAME = "results_bernoulli"

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


    const.BERNOULLI_FEED0_CPS_MTX = lambda ref, RE: const.ASSEMBLER_COMPRESS_CTRL_NAME_MTX(ref, "Kinf", RE)

    build_mesh(const, refs)
    solve_newton(const, refs, REs)
    REs = REscompute
    assemble_lqr(const, refs, REs)
    simulate(const, refs, REs)
    solve_bernoulli(const, refs, REs)
    # control(const, refs, REs)

    # control Makes only sense if bernoulli feedback is available
    for ref in refs:
        for RE in REs:
            try:
                print "Control ref={0:d} RE={1:d}".format(ref, RE)
                linearizedctrl = LinearizedCtrl(const, ref, RE)
                linearizedctrl.solve_ode()
                linearizedctrl.save_log()
            except:
                print "An exception in control ref={0:d} RE={1:d}".format(ref, RE)
                print traceback.format_exc()
                continue

    # compute_eigen(const, refs, [max(REs)])
    # plot(const, refs, REs)
