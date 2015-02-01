import src.karman_const as const
const.OUTPUTDIR_NAME = "results_default"
const.ASSEMBLER_OBSERVER_POINT1_X = 2.5
const.ASSEMBLER_OBSERVER_POINT1_Y = 0.5
const.ASSEMBLER_OBSERVER_POINT2_X = 3.0
const.ASSEMBLER_OBSERVER_POINT2_Y = 0.5

from src.mesh.karman import MeshBuilder
from src.stationary.newton import Newton
from src.lqr.assembler import Assembler
from src.lqr.compress_assembler import CompressAssembler
from src.lqr.bernoulli import Bernoulli
from src.lqr.lqr_solver import LQR_Solver
from src.lqr.linearized_sim import LinearizedSim
from src.lqr.linearized_ctrl import LinearizedCtrl


def build_mesh(ref):
    print "build mesh"
    meshbuilder = MeshBuilder()
    for r in range(0, ref):
        meshbuilder.refine()
    meshbuilder.save()


def solve_newton(ref, REs):
    print "solve newton"
    REinitial = None
    for RE in REs:
        newton = Newton(ref, RE, REinitial)
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
    ref = 3
    REs = range(100, 700, 100)
    dt = 0.005
    T = 30
    pertubationeps = 0.25

    build_mesh(ref)
    solve_newton(ref, REs)
    lqr_assembler(ref, REs)
    bernoulli(ref, REs)
    lqr_solver(ref, REs)
    sim_and_control(ref, REs, pertubationeps, dt, T)