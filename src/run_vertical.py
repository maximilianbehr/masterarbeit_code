import src.karman_const as const
from src.mesh.karman import MeshBuilder
from src.stationary.newton import Newton
from src.lqr.assembler import Assembler
from src.lqr.compress_assembler import CompressAssembler
from src.lqr.bernoulli import Bernoulli
from src.lqr.lqr_solver import LQR_Solver
from src.lqr.linearized_sim import LinearizedSim
from src.lqr.linearized_ctrl import LinearizedCtrl
from src.lqr.plotter import Plotter


def build_mesh(ref):
    print "build mesh"
    meshbuilder = MeshBuilder()
    for r in range(0, ref):
        meshbuilder.refine()
    meshbuilder.save()


def solve_newton(const, ref, REs):
    print "solve newton"
    REinitial = None
    for RE in REs:
        newton = Newton(const, ref, RE, REinitial)
        newton.solve()
        newton.save()
        REinitial = RE


def lqr_assembler(ref, REs):
    print "lqr assembler"
    for RE in REs:
        assembler = Assembler(ref, RE)
        assembler.save_mat()
        assembler.save_mtx()
        simassembler = CompressAssembler(ref, RE, "sim")
        simassembler.save()
        ctrlassembler = CompressAssembler(ref, RE, "ctrl")
        ctrlassembler.save()


def bernoulli(ref, REs):
    print "bernoulli"
    for RE in REs:
        bernoulli = Bernoulli(ref, RE)
        bernoulli.solve()
        bernoulli.save()


def lqr_solver(ref, REs):
    print "lqr solver"
    for RE in REs:
        lqrsolver = LQR_Solver(ref, RE)
        lqrsolver.solve()
        lqrsolver.save()


def sim_and_control(ref, REs, pertubationeps, dt, T):
    print "sim and control"
    for RE in REs:
        linearized = LinearizedSim(ref, RE, pertubationeps, dt, T)
        linearized.solve_ode()
        linearized.save_log()
        linearized = LinearizedCtrl(ref, RE, pertubationeps, dt, T)
        linearized.solve_ode()
        linearized.save_log()

if __name__ == "__main__":

    ref = 4
    REs = range(100, 1500, 100)
    dt = 0.001
    T = 60
    pertubationeps = 0.25

    for pre in [0, 1, 2, 3, 4, 5]:
        const.OUTPUTDIR_NAME = "results_vertical_pre_{0:d}s".format(pre)
        const.LINEARIZED_CTRL_START_CONTROLLING = pre

        build_mesh(ref)
        solve_newton(const, ref, REs)
        lqr_assembler(ref, REs)
        bernoulli(ref, REs)
        lqr_solver(ref, REs)
        sim_and_control(ref, REs, pertubationeps, dt, T)
        plotter = Plotter(ref, REs)
        plotter.plot_linearized_sim1()
        plotter.plot_linearized_ctrl1()
        plotter.plot_linearized_ctrl2()
        plotter.plot_linearized_ctrl3()
