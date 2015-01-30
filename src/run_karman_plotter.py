from src.lqr.plotter import Plotter
import traceback


REs = range(100, 600, 100)
refs = [3]

for ref in refs:
        try:
            plotter = Plotter(ref, REs)
            plotter.plot_linearized_sim()

        except Exception, e:
            print traceback.print_exc()


