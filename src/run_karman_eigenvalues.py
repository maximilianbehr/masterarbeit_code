import src.karman_const as const
from src.lqr.eigen import Eigen
import traceback
from src.aux import gettime

# set Reynoldsnumbers and refinements
REs = range(100, 500, 100)
refs = [3]


for ref in refs:
    for RE in REs:
        print "ref={0:d} RE={1:d}".format(ref, RE)
        eig = Eigen(const, ref, RE)

        print "{0:s}, start Eigenvalue System".format(gettime())
        eig.compute_eig_sys()

        print "{0:s}, start Eigenvalue Riccati".format(gettime())
        eig.compute_eig_ric()

        print "{0:s}, start Eigenvalue Bernoulli".format(gettime())
        eig.compute_eig_ber()
