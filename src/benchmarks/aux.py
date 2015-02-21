from src.stationary.newton import Newton
from src.lqr.assembler import Assembler
from src.lqr.compress_assembler import CompressAssembler
from src.lqr.bernoulli import Bernoulli
from src.lqr.lqr_solver import LQR_Solver
from src.lqr.linearized_sim import LinearizedSim
from src.lqr.linearized_ctrl import LinearizedCtrl
from src.lqr.plotter import Plotter
from src.lqr.eigen import Eigen
from src.aux import *

def build_mesh(const, refs):
    if const.NAME == "bws":
        from src.benchmarks.bws.mesh.bws import MeshBuilder
    elif const.NAME == "karman":
        from src.benchmarks.karman.mesh.karman import MeshBuilder

    # generate a mesh builder
    meshbuilder = MeshBuilder(const)

    for r in range(0, max(refs)):
        print "Build Mesh ref={0:d}".format(r)
        # refine and save mesh
        meshbuilder.refine()
        meshbuilder.save()


def solve_newton(const, refs, REs):
    for ref in refs:
        REinitial = None
        for RE in REs:
            print "Newton ref={0:d} RE={1:d}".format(ref, RE)
            newton = Newton(const, ref, RE, REinitial)
            newton.solve()
            newton.save()
            REinitial = RE


def assemble_lqr(const, refs, REs):
    for ref in refs:
        for RE in REs:
            print "Assemble ref={0:d} RE={1:d}".format(ref, RE)
            assembler = Assembler(const, ref, RE)
            assembler.save_mat()
            assembler.save_mtx()
            simassembler = CompressAssembler(const, ref, RE, "sim")
            simassembler.save_mtx()
            ctrlassembler = CompressAssembler(const, ref, RE, "ctrl")
            ctrlassembler.save_mtx()


def solve_bernoulli(const, refs, REs):
    for ref in refs:
        for RE in REs:
            print "Bernoulli ref={0:d} RE={1:d}".format(ref, RE)
            bernoulli = Bernoulli(const, ref, RE)
            bernoulli.solve()
            bernoulli.save()


def compute_eigen(const, refs, REs):
    for ref in refs:
        for RE in REs:
            print "Eigen ref={0:d} RE={1:d}".format(ref, RE)
            clear_prof_data()
            eig = Eigen(const, ref, RE)

            eig.compute_eig_sys()
            print_prof_data()
            clear_prof_data()

            eig.compute_eig_ric()
            print_prof_data()
            clear_prof_data()

            eig.compute_eig_ber()
            print_prof_data()
            clear_prof_data()

            eig.save()
            eig.plot()


def solve_lqr(const, refs, REs):
    for ref in refs:
        for RE in REs:
            print "LQR Solver ref={0:d} RE={1:d}".format(ref, RE)
            lqrsolver = LQR_Solver(const, ref, RE)
            lqrsolver.solve()
            lqrsolver.save()


def simulate(const, refs, REs):
    for ref in refs:
        for RE in REs:
            print "Simulate ref={0:d} RE={1:d}".format(ref, RE)
            linearizedsim = LinearizedSim(const, ref, RE)
            linearizedsim.solve_ode()
            linearizedsim.save_log()


def control(const, refs, REs):
    for ref in refs:
        for RE in REs:
            print "Control ref={0:d} RE={1:d}".format(ref, RE)
            linearizedctrl = LinearizedCtrl(const, ref, RE)
            linearizedctrl.solve_ode()
            linearizedctrl.save_log()

def plot(const, refs, REs):
    for ref in refs:
        print "Plotter ref ={0:d}".format(ref)
        plotter = Plotter(const, ref, REs)
        plotter.plot_linearized_sim()
        plotter.plot_linearized_ctrl()
        plotter.plot_linearized_ctrl_controls()
        plotter.plot_linearized_ctrl_output()
