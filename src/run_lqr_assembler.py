from src.lqr.assembler import Assembler
import traceback

# karman
REs = range(100, 5000, 100)
refs = [1, 2, 3]

for ref in refs:
    for RE in REs:

        try:
            print "ref={0:d} RE={1:d}".format(ref,RE)
            assembler = Assembler(ref, RE)
            assembler.lns_variational()
            assembler.lns_ublas()
            assembler.lns_npsc()
            print "save mtx"
            assembler.save_lns_mtx()
            assembler.save_lns_mat()
        except Exception, e:
            traceback.print_exc()



