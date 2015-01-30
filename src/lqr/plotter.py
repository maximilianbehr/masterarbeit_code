import src.karman_const as const
import numpy as np
import matplotlib.pyplot as plt
import os



class Plotter():

    def __init__(self, ref, REs):

        # set parameters
        self.ref = ref
        self.REs = REs

    def plot_linearized_sim(self):

        for RE in self.REs:
            log = np.loadtxt(const.LINEARIZED_SIM_LOG(self.ref, RE))
            plt.semilogy(log[:, 0], log[:, 1], label="$Re={0:d}$".format(RE))

        plt.title("Simulation")
        plt.grid(True)

        plt.xlabel("$t \mathrm{\;in\;Sekunden}$")
        plt.ylabel("$\|\| u_{\delta} \|\|_2$")
        plt.xticks(np.arange(min(log[:, 0]), max(log[:, 0])+1, 1))

        # legend for each plot
        plt.legend(bbox_to_anchor=(0.0, 0.0), loc=3)

        #plt.show(block=True)
        if not os.path.exists(os.path.dirname(const.PLOTTER_LINEARIZED_SIM_LOG_EPS(self.ref))):
                os.makedirs(os.path.dirname(const.PLOTTER_LINEARIZED_SIM_LOG_EPS(self.ref)))

        plt.savefig(const.PLOTTER_LINEARIZED_SIM_LOG_EPS(self.ref))
        plt.savefig(const.PLOTTER_LINEARIZED_SIM_LOG_PNG(self.ref))
        plt.savefig(const.PLOTTER_LINEARIZED_SIM_LOG_JPEG(self.ref))






