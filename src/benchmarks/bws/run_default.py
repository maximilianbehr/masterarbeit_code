# -*- coding: utf-8 -*-
import src.benchmarks.bws.bws_const as const
from src.benchmarks.aux import *

if __name__ == "__main__":
    # set Reynoldsnumbers and refinements and Parameters
    REs = range(100, 1000, 100)
    # refs = [1, 2, 3, 4, 5]
    #refs = [1]
    refs = [2]

    const.LINEARIZED_SIM_T = 60.0
    const.LINEARIZED_SIM_DT = 0.0025
    const.LINEARIZED_CTRL_T = 60.0
    const.LINEARIZED_CTRL_DT = 0.0025
    const.LINEARIZED_SIM_INFO = 0.1
    const.LINEARIZED_CTRL_INFO = 0.1
    const.LQR_ADI_OUTPUT = 1

    # build mesh
    for lam in [-0.5, 0.0, 0.5]:


        const.STATIONARY_CONTROL_UPPER.lam = lam
        const.STATIONARY_CONTROL_LOWER.lam = lam
        const.OUTPUTDIR_NAME = "results_bws_lam_{0:3f}".format(lam)

        # build mesh
        build_mesh(const, refs)

        # solve stationary
        solve_newton(const, refs, REs)

        # assemble lqr
        REs = [200,900]
        assemble_lqr(const, refs, REs)

        # solve bernoulli
        solve_bernoulli(const, refs, REs)

        # solve lqr
        solve_lqr(const, refs, REs)

        # simulate
        simulate(const, refs, REs)

        # simulate with control
        control(const, refs, REs)

        # compute eigenvalues
        #REs = [max(REs)]
        #compute_eigen(const, refs, REs)


