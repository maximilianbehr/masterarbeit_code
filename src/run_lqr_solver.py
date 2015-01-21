from src.lqr.lqr_solver import LQR_Solver
import traceback

REs = range(100, 5000, 100)
refs = [2]

for ref in refs:
    REinitial = None
    for RE in REs:

        print "ref = {0:d} RE = {1:d}".format(ref, RE)
        try:
            lqrsolver = LQR_Solver(ref, RE, REinitial)
            lqrsolver.solve()
            print "save"
            lqrsolver.save()
            REinitial = RE
        except Exception, e:
            print traceback.print_exc()
            # reynoldsnumber to high break
            break
