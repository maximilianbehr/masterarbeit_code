
from src.lqr.linearized_ctrl import LinearizedCtrl
import traceback


REs = range(900, 1000, 100)
refs = [3]
dt = 0.01
T = 30
pertubationeps = 0.25

for ref in refs:
    for RE in REs:
        print "ref={0:d} RE={1:d}".format(ref, RE)
        try:
            linearized = LinearizedCtrl(ref, RE, pertubationeps, dt, T)
            linearized.build()
            linearized.solve_ode()
        except Exception, e:
            print traceback.print_exc()


