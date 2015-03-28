from src.amg.solver_diagnostics import solver_diagnostics
import scipy.io as scio
import scipy.sparse as scsp
import src.benchmarks.karman.karman_const as const

#parameters
ref = 1
RE = 20
dt = 0.001
const.OUTPUTDIR_NAME = "results_karman_ref_1"

# read compress system for simuation
names = ["M", "M_BOUNDARY_CTRL", "S", "R", "K", "G", "GT", "B", "C"]
mat = {}
for name in names:
    mat[name] = scio.mmread(const.ASSEMBLER_COMPRESS_SIM_NAME_MTX(ref, name, RE))


nv, np = mat["G"].shape


# build system matrices
Asys = -mat["S"]-mat["R"]-mat["K"]
u = scsp.hstack([mat["M"] - dt*Asys, dt*(-mat["G"])])
l = scsp.hstack([dt * (-mat["GT"]), scsp.csr_matrix((np, np))])
Msys_ode = scsp.vstack([u, l]).tocsc()

solver_diagnostics(Msys_ode.tocsr(), fname='sim_karman_{0:d}'.format(RE),\
    cycle_list=['V', 'W'], \
    definiteness='indefinite', \
    symmetry='nonsymmetric', coarse_size_list = [ (300, 'pinv'), (300, 'splu') ])






