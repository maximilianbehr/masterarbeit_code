from pycmess import options
from pycmess import equation_dae2
from pycmess import PYCMESS_OP_TRANSPOSE
from pycmess import lrnm
from pycmess import lrnm_res

# create opt instance
opt = options()
opt.adi.output = 1;
opt.nm.output = 1;
opt.nm.res2_tol = 1e-5;

eqn = equation_dae2()
eqn.M = mat_numpyscipy["M"]
eqn.A = 1.0 / nu * (
    -mat_numpyscipy["S"] - 1.0 / penalty_eps * (mat_numpyscipy["M_lower"] + mat_numpyscipy["M_upper"])) - \
        mat_numpyscipy[
            "K"] - mat_numpyscipy["R"]
eqn.G = -1 * mat_numpyscipy["G"]
eqn.B = 1.0 / nu * 1.0 / penalty_eps * mat_numpyscipy["B"]
eqn.C = Ccsr
eqn.delta = -0.02

# solve generalized riccati equation
opt.adi.type = PYCMESS_OP_TRANSPOSE
result = lrnm(eqn, opt)
Z = result[0]
res2 = result[1]
iter = result[2]
print res2[res2.size - 1]
result = lrnm_res(eqn, opt, Z)
res = result[0]
rel = result[1]
print "rel= %e \t res= %e \t res/rel= %e \n" % (rel, res, res / rel)