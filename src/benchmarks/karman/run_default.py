# -*- coding: utf-8 -*-
import src.benchmarks.karman.karman_const as const
from src.benchmarks.aux import build_mesh
from src.benchmarks.aux import solve_newton
from src.benchmarks.aux import assemble_lqr
from src.benchmarks.aux import solve_bernoulli
from src.benchmarks.aux import simulate

if __name__ == "__main__":
    # set Reynoldsnumbers and refinements and Parameters
    REs = range(100, 1000, 100)
    refs = [3]

    const.LINEARIZED_SIM_DT = 0.005
    const.LINEARIZED_SIM_T = 30

    # build mesh
    build_mesh(const, refs)

    # solve stationary
    solve_newton(const, refs, REs)

    # assemble lqr
    assemble_lqr(const, refs, REs)

    # simulate
    simulate(const, refs, REs)

    # solve bernoulli
    # solve_bernoulli(const, refs, REs)
