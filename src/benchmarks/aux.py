from src.stationary.newton import Newton
from src.lqr.assembler import Assembler
from src.lqr.compress_assembler import CompressAssembler
from src.lqr.bernoulli import Bernoulli
from src.lqr.lqr_solver import LQR_Solver
from src.lqr.linearized_sim_petsc import LinearizedSimPETSC
from src.lqr.linearized_ctrl import LinearizedCtrl
from src.lqr.plotter import Plotter
from src.lqr.eigen import Eigen
from src.aux import *
import traceback

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
            if REinitial:
                winitial = const.STATIONARY_W_XML(ref, REinitial)
                newton = Newton(const, ref, RE, winitial)
            else:
                newton = Newton(const, ref, RE)
            try:
                newton.solve()
                newton.save()
                REinitial = RE
            except:
                print "An exception in solve newton ref={0:d} RE={1:d}".format(ref, RE)
                print traceback.format_exc()
                break

def assemble_lqr(const, refs, REs):
    for ref in refs:
        for RE in REs:
            try:
                print "Assemble ref={0:d} RE={1:d}".format(ref, RE)
                assembler = Assembler(const, ref, RE)
                # assembler.save_mat()
                assembler.save_mtx()
                print "SimAssemble ref={0:d} RE={1:d}".format(ref, RE)
                simassembler = CompressAssembler(const, ref, RE, "sim")
                simassembler.save_mtx()
                print "CtrlAssemble ref={0:d} RE={1:d}".format(ref, RE)
                ctrlassembler = CompressAssembler(const, ref, RE, "ctrl")
                ctrlassembler.save_mtx()
            except:
                print "An exception in assemble_lqr ref={0:d} RE={1:d}".format(ref, RE)
                print traceback.format_exc()
                break


def solve_bernoulli(const, refs, REs):
    for ref in refs:
        for RE in REs:
            try:
                print "Bernoulli ref={0:d} RE={1:d}".format(ref, RE)
                bernoulli = Bernoulli(const, ref, RE)
                bernoulli.solve()
                bernoulli.save()
            except:
                print "An exception in solve_bernoulli ref={0:d} RE={1:d}".format(ref, RE)
                print traceback.format_exc()
                break

def compute_eigen(const, refs, REs):
    for ref in refs:
        for RE in REs:
            try:
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
            except:
                print "An exception in compute eigen ref={0:d} RE={1:d}".format(ref, RE)
                print traceback.format_exc()
                break

def solve_lqr(const, refs, REs):
    for ref in refs:
        for RE in REs:
            try:
                print "LQR Solver ref={0:d} RE={1:d}".format(ref, RE)
                lqrsolver = LQR_Solver(const, ref, RE)
                lqrsolver.solve()
                lqrsolver.save()
            except:
                print "An exception in solve_lqr ref={0:d} RE={1:d}".format(ref, RE)
                print traceback.format_exc()
                break

def simulate(const, refs, REs):
    for ref in refs:
        for RE in REs:
            try:
                print "Simulate ref={0:d} RE={1:d}".format(ref, RE)
                # linearizedsim = LinearizedSim(const, ref, RE)
                linearizedsim = LinearizedSimPETSC(const, ref, RE)
                linearizedsim.solve_ode()
                linearizedsim.save_log()
            except:
                print "An exception in simulate ref={0:d} RE={1:d}".format(ref, RE)
                print traceback.format_exc()
                break

def control(const, refs, REs):
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
                break


def plot(const, refs, REs):
    for ref in refs:
        try:
            print "Plotter ref ={0:d}".format(ref)
            plotter = Plotter(const, ref, REs)
            plotter.plot_linearized_sim()
            plotter.plot_linearized_ctrl()
            plotter.plot_linearized_ctrl_controls()
            plotter.plot_linearized_ctrl_output()
        except:
            print "An exception in plot ref={0:d}".format(ref)
            print traceback.format_exc()
            break
