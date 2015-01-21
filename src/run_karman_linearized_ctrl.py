
from src.lqr.linearized_ctrl import LinearizedCtrl
import traceback


REs = range(100, 1000, 100)
refs = [2]
dt = 0.005
T = 30
pertubationeps = 0.5

for ref in refs:
    for RE in REs:
        print "ref={0:d} RE={1:d}".format(ref, RE)
        try:
            linearized = LinearizedCtrl(ref, RE, pertubationeps, dt, T)
            linearized.build()
            linearized.solve_ode()
        except Exception, e:
            print traceback.print_exc()


