import numpy as np
import matplotlib.pyplot as plt
from src.aux import createdir
import matplotlib


class Plotter():

    def __init__(self, const, ref, REs):

        # set parameters
        self.ref = ref
        self.REs = REs
        self.const = const
        self.legendsize = 6

        # set random colors
        hexcolors = matplotlib.colors.cnames.values()
        idhex = np.linspace(10, len(hexcolors)-10, len(REs))
        hexlist = idhex.astype(int).tolist()
        # self.color = [(rand(), rand(), rand()) for i in range(1, 20)]
        self.color = [hexcolors[hex] for hex in hexlist]


    def plot_linearized_sim(self):
        # full plot
        self._plot_linearized_sim(0.0, -1)
        percents = np.linspace(0.0, 0.6, 3)
        for num in range(percents.shape[0]):
            self._plot_linearized_sim(percents[num], num)

    def plot_linearized_ctrl(self):
        # full plot
        self._plot_linearized_ctrl(0.0, -1)
        percents = np.linspace(0.0, 0.6, 3)
        for num in range(percents.shape[0]):
            self._plot_linearized_ctrl(percents[num], num)

    def plot_linearized_ctrl_controls(self):
        # full plot
        self._plot_linearized_ctrl_controls(0.0, -1)
        percents = np.linspace(0.0, 0.6, 3)
        for num in range(percents.shape[0]):
            self._plot_linearized_ctrl_controls(percents[num], num)

    def plot_linearized_ctrl_output(self):
        # full plot
        self._plot_linearized_ctrl_controls(0.0, -1)
        percents = np.linspace(0.0, 0.6, 3)
        for num in range(percents.shape[0]):
            self._plot_linearized_ctrl_output(percents[num], num)


    def _plot_linearized_sim(self, percent, num):

        k = 0
        for RE in self.REs:
            log = np.loadtxt(self.const.LINEARIZED_SIM_LOG(self.ref, RE))
            length = log.shape[0]
            plt.semilogy(log[int(percent*length):, 0], log[int(percent*length):, 1],\
                         label="$Re={0:d}$".format(RE), color=self.color[k])
            k += 1

        plt.title("Simulation $v_{\delta}$")
        plt.grid(True)
        plt.xlabel("$t \mathrm{\;in\;Sekunden}$")
        plt.ylabel("$\|\| v_{\delta} \|\|_2$")
        l = int(min(log[int(percent*length):, 0]))
        u = int(max(log[int(percent*length):, 0]))+1
        plt.xticks(np.arange(l, u, 4))
        # legend for each plot
        plt.legend(bbox_to_anchor=(0.0, 0.0), loc=3, prop={'size': self.legendsize})
        #plt.show(block=True)
        createdir(self.const.PLOTTER_LINEARIZED_SIM_LOG(self.ref, num, "png"))
        plt.savefig(self.const.PLOTTER_LINEARIZED_SIM_LOG(self.ref, num, "png"))
        plt.savefig(self.const.PLOTTER_LINEARIZED_SIM_LOG(self.ref, num, "jpeg"))
        plt.savefig(self.const.PLOTTER_LINEARIZED_SIM_LOG(self.ref, num, "eps"))
        plt.close()

    def _plot_linearized_ctrl(self, percent, num):

        k = 0
        for RE in self.REs:
            log = np.loadtxt(self.const.LINEARIZED_CTRL_LOG(self.ref, RE))
            length = log.shape[0]
            plt.semilogy(log[int(percent*length):, 0], log[int(percent*length):, 1],\
                         label="$Re={0:d}$".format(RE), color=self.color[k])
            k += 1

        plt.title("Simulation mit Steuerung $v_{\delta}$")
        plt.grid(True)
        plt.xlabel("$t \mathrm{\;in\;Sekunden}$")
        plt.ylabel("$\|\| v_{\delta} \|\|_2$")
        l = int(min(log[int(percent*length):, 0]))
        u = int(max(log[int(percent*length):, 0]))+1
        plt.xticks(np.arange(l, u, 4))
        # legend for each plot
        plt.legend(bbox_to_anchor=(0.0, 0.0), loc=3, prop={'size': self.legendsize})
        #plt.show(block=True)
        createdir(self.const.PLOTTER_LINEARIZED_CTRL_LOG(self.ref, num, "png"))
        plt.savefig(self.const.PLOTTER_LINEARIZED_CTRL_LOG(self.ref, num, "png"))
        plt.savefig(self.const.PLOTTER_LINEARIZED_CTRL_LOG(self.ref, num, "jpeg"))
        plt.savefig(self.const.PLOTTER_LINEARIZED_CTRL_LOG(self.ref, num, "eps"))
        plt.close()

    def _plot_linearized_ctrl_controls(self, percent, num):
        k = 0
        for RE in self.REs:
            log = np.loadtxt(self.const.LINEARIZED_CTRL_LOG(self.ref, RE))
            length = log.shape[0]
            plt.plot(log[int(percent*length):, 0], log[int(percent*length):, 2],\
                         "-", color=self.color[k])
            plt.plot(log[int(percent*length):, 0], log[int(percent*length):, 3],\
                         "--", label="$Re={0:d}$".format(RE), color=self.color[k])
            k += 1

        plt.title("Simulation mit Steuerung ")
        plt.grid(True)
        plt.xlabel("$t \mathrm{\;in\;Sekunden}$")
        plt.ylabel("$Cv_{\delta}$")
        l = int(min(log[int(percent*length):, 0]))
        u = int(max(log[int(percent*length):, 0]))+1
        plt.xticks(np.arange(l, u, 4))
        # legend for each plot
        plt.legend(bbox_to_anchor=(0.0, 0.0), loc=3, prop={'size': self.legendsize})
        #plt.show(block=True)
        createdir(self.const.PLOTTER_LINEARIZED_CTRL_CONTROLS_LOG(self.ref, num, "png"))
        plt.savefig(self.const.PLOTTER_LINEARIZED_CTRL_CONTROLS_LOG(self.ref, num, "png"))
        plt.savefig(self.const.PLOTTER_LINEARIZED_CTRL_CONTROLS_LOG(self.ref, num, "jpeg"))
        plt.savefig(self.const.PLOTTER_LINEARIZED_CTRL_CONTROLS_LOG(self.ref, num, "eps"))
        plt.close()

    def _plot_linearized_ctrl_output(self, percent, num):
        k = 0
        for RE in self.REs:
            log = np.loadtxt(self.const.LINEARIZED_CTRL_LOG(self.ref, RE))
            length = log.shape[0]
            plt.plot(log[int(percent*length):, 0], log[int(percent*length):, 4],\
                         "-", color=self.color[k])
            plt.plot(log[int(percent*length):, 0], log[int(percent*length):, 5],\
                         "--", label="$Re={0:d}$".format(RE), color=self.color[k])
            k += 1

        plt.title("Simulation mit Steuerung $v_delta$")
        plt.grid(True)
        plt.xlabel("$t \mathrm{\;in\;Sekunden}$")
        plt.ylabel("$u_1,u_2$")
        l = int(min(log[int(percent*length):, 0]))
        u = int(max(log[int(percent*length):, 0]))+1
        plt.xticks(np.arange(l, u, 4))
        # legend for each plot
        plt.legend(bbox_to_anchor=(0.0, 0.0), loc=3, prop={'size': self.legendsize})
        #plt.show(block=True)
        createdir(self.const.PLOTTER_LINEARIZED_CTRL_OUTPUT_LOG(self.ref, num, "png"))
        plt.savefig(self.const.PLOTTER_LINEARIZED_CTRL_OUTPUT_LOG(self.ref, num, "png"))
        plt.savefig(self.const.PLOTTER_LINEARIZED_CTRL_OUTPUT_LOG(self.ref, num, "jpeg"))
        plt.savefig(self.const.PLOTTER_LINEARIZED_CTRL_OUTPUT_LOG(self.ref, num, "eps"))
        plt.close()

