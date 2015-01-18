from src.lqr.assembler import Assembler


# karman
REs = range(100, 3000, 100)
refs = [2, 3, 4]

for ref in refs:
    for RE in REs:

        assembler = Assembler(ref, RE)
        assembler.lns_variational()
        assembler.lns_ublas()
        assembler.save_lns_mtx()
        assembler.save_lns_mat()



