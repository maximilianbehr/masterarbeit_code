from src.lqr.assembler import Assembler
from src.lqr.compress_assembler import CompressAssembler
import traceback

# karman
REs = range(100, 1000, 100)
refs = [3]

for ref in refs:
    for RE in REs:

        try:
            print "ref={0:d} RE={1:d}".format(ref,RE)
            assembler = Assembler(ref, RE)
            assembler.lns_variational()
            assembler.lns_ublas()
            assembler.lns_npsc()

            print "save lqr mtx"
            assembler.save_lns_mtx()
            assembler.save_lns_mat()

            print "save compress sim mtx"
            simassembler = CompressAssembler(ref, RE, "sim")
            simassembler.save()

            print "save compress ctrl mtx"
            ctrlassembler = CompressAssembler(ref, RE, "ctrl")
            ctrlassembler.save()

        except Exception, e:
            traceback.print_exc()



