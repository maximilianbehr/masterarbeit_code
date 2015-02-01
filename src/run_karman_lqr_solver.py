from src.lqr.lqr_solver import LQR_Solver
import traceback

REs = range(100, 700, 100)
refs = [3]

for ref in refs:
    for RE in REs:

        print "ref = {0:d} RE = {1:d}".format(ref, RE)
        try:
            lqrsolver = LQR_Solver(ref, RE)
            lqrsolver.solve()
            lqrsolver.save()
        except Exception, e:
            print traceback.print_exc()
            # reynoldsnumber to high break
            #break
