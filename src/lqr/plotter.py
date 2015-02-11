import src.karman_const as const
import numpy as np
import matplotlib.pyplot as plt
import os
from src.aux import createdir
from numpy.random import rand


class Plotter():

    def __init__(self, ref, REs):

        # set parameters
        self.ref = ref
        self.REs = REs
        # set random colors
        self.color = [(rand(), rand(), rand()) for i in range(1, 20)]

    def plot_linearized_sim1(self):

        k = 0
        for RE in self.REs:
            log = np.loadtxt(const.LINEARIZED_SIM_LOG(self.ref, RE))
            plt.semilogy(log[:, 0], log[:, 1], label="$Re={0:d}$".format(RE), color=self.color[k])
            k += 1

        plt.title("Simulation")
        plt.grid(True)

        plt.xlabel("$t \mathrm{\;in\;Sekunden}$")
        plt.ylabel("$\|\| u_{\delta} \|\|_2$")
        plt.xticks(np.arange(min(log[:, 0]), max(log[:, 0])+1, 1))

        # legend for each plot
        plt.legend(bbox_to_anchor=(0.0, 0.0), loc=3)

        #plt.show(block=True)
        createdir(const.PLOTTER_LINEARIZED_SIM_LOG1(self.ref, "png"))

        plt.savefig(const.PLOTTER_LINEARIZED_SIM_LOG1(self.ref, "png"))
        plt.savefig(const.PLOTTER_LINEARIZED_SIM_LOG1(self.ref, "jpeg"))
        plt.savefig(const.PLOTTER_LINEARIZED_SIM_LOG1(self.ref, "eps"))

        plt.close()

    def plot_linearized_ctrl1(self):

        k = 0
        for RE in self.REs:
            log = np.loadtxt(const.LINEARIZED_CTRL_LOG(self.ref, RE))
            plt.semilogy(log[:, 0], log[:, 1], label="$Re={0:d}$".format(RE), color=self.color[k])
            k += 1

        plt.title("Simulation")
        plt.grid(True)

        plt.xlabel("$t \mathrm{\;in\;Sekunden}$")
        plt.ylabel("$\|\| u_{\delta} \|\|_2$")
        plt.xticks(np.arange(min(log[:, 0]), max(log[:, 0])+1, 2))

        # legend for each plot
        plt.legend(bbox_to_anchor=(0.0, 0.0), loc=3)

        #plt.show(block=True)
        createdir(const.PLOTTER_LINEARIZED_CTRL_LOG1(self.ref, "png"))
        plt.savefig(const.PLOTTER_LINEARIZED_CTRL_LOG1(self.ref, "png"))
        plt.savefig(const.PLOTTER_LINEARIZED_CTRL_LOG1(self.ref, "jpeg"))
        plt.savefig(const.PLOTTER_LINEARIZED_CTRL_LOG1(self.ref, "eps"))
        plt.close()

    def plot_linearized_ctrl2(self):

        k = 0
        for RE in self.REs:
            log = np.loadtxt(const.LINEARIZED_CTRL_LOG(self.ref, RE))
            plt.semilogy(log[:, 0], log[:, 2], "-", color=self.color[k])
            plt.semilogy(log[:, 0], log[:, 3], "--", label="$Re={0:d}$".format(RE), color=self.color[k])
            k += 1

        plt.title("Simulation")
        plt.grid(True)

        plt.xlabel("$t \mathrm{\;in\;Sekunden}$")
        plt.ylabel("$u_{c1}, u_{c2}$")
        plt.xticks(np.arange(min(log[:, 0]), max(log[:, 0])+1, 2))

        # legend for each plot
        plt.legend(bbox_to_anchor=(0.0, 0.0), loc=3)

        #plt.show(block=True)
        createdir(const.PLOTTER_LINEARIZED_CTRL_LOG2(self.ref, "png"))
        plt.savefig(const.PLOTTER_LINEARIZED_CTRL_LOG2(self.ref, "png"))
        plt.savefig(const.PLOTTER_LINEARIZED_CTRL_LOG2(self.ref, "jpeg"))
        plt.savefig(const.PLOTTER_LINEARIZED_CTRL_LOG2(self.ref, "eps"))
        plt.close()


    def plot_linearized_ctrl3(self):
        k = 0
        for RE in self.REs:
            log = np.loadtxt(const.LINEARIZED_CTRL_LOG(self.ref, RE))
            plt.plot(log[:, 0], log[:, 2], "-", color=self.color[k])
            plt.plot(log[:, 0], log[:, 3], "--", label="$Re={0:d}$".format(RE))
            k += 1

        plt.title("Simulation")
        plt.grid(True)

        plt.xlabel("$t \mathrm{\;in\;Sekunden}$")
        plt.ylabel("$u_{c1}, u_{c2}$")
        plt.xticks(np.arange(min(log[:, 0]), max(log[:, 0])+1, 2))

        # legend for each plot
        plt.legend(bbox_to_anchor=(0.0, 0.0), loc=3)

        #plt.show(block=True)
        createdir(const.PLOTTER_LINEARIZED_CTRL_LOG3(self.ref, "png"))
        plt.savefig(const.PLOTTER_LINEARIZED_CTRL_LOG3(self.ref, "png"))
        plt.savefig(const.PLOTTER_LINEARIZED_CTRL_LOG3(self.ref, "jpeg"))
        plt.savefig(const.PLOTTER_LINEARIZED_CTRL_LOG3(self.ref, "eps"))
        plt.close()

"""
        #split set of eigenvalues in stable, unstable and zeros (hopefully no)
        stable_eigs = eigs[eigs.real<0]
        unstable_eigs = eigs[eigs.real>0]
        zero_eigs = eigs[eigs.real==0]

        #plot eigenvalues
        fig, ax = plt.subplots()
        ax.plot(stable_eigs.real, stable_eigs.imag, "rx")
        ax.plot(unstable_eigs.real,unstable_eigs.imag, "bx")
        ax.plot(zero_eigs.real,zero_eigs.imag,"gx")
        plt.axvline(x=1.0/self.delta, linewidth=1, color="g",ls="dashed")
        xlimit = numpy.max(numpy.ceil(numpy.absolute(eigs.real)))
        ylimit = numpy.max(numpy.ceil(numpy.absolute(eigs.imag)))
        plt.xlim((-xlimit, xlimit))
        plt.ylim((-ylimit, ylimit))
        plt.xscale("symlog")
        plt.xlabel("Real")
        plt.ylabel("Imaginary")
        #plt.show()
        plt.savefig(self.options["eig_eps"])
        plt.close("all")
"""