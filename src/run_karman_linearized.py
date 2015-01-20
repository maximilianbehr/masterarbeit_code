
from src.lqr.linearized import Linearized


REs = range(100, 1500, 200)
#REs = [200, 400, 500]
refs = [3]
dt = 0.01
T = 20
pertubationeps = 0.5

for ref in refs:
    for RE in REs:
        print "ref={0:d} RE={1:d}".format(ref, RE)
        linearized = Linearized(ref, RE, pertubationeps, dt, T)
        linearized.solve_ode()



