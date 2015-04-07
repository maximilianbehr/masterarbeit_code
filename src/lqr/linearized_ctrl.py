import scipy.io as scio
import scipy.sparse as scsp
import scipy.sparse.linalg as scspli
import scipy.linalg as scla
from dolfin import *
import numpy as np
from src.aux import gettime
from src.lqr.linearized_aux import smw, stable_timestep

# set floating point error handling
np.seterr(all="raise", divide="raise", over="raise", under="raise", invalid="raise")


class LinearizedCtrl():
    def __init__(self, const, ref, RE):

        # parameters
        self.ref = ref
        self.RE = RE
        self.const = const
        self.T = self.const.LINEARIZED_CTRL_T
        self.pertubationeps = self.const.LINEARIZED_CTRL_PERTUBATIONEPS

        # mesh and function spaces
        self.mesh = Mesh(const.MESH_XML(ref))
        self.V = VectorFunctionSpace(self.mesh, const.V, const.V_DIM)
        self.Q = FunctionSpace(self.mesh, const.Q, const.Q_DIM)
        self.u_stat = Function(self.V, const.STATIONARY_U_XML(ref, RE))

        if const.LINEARIZED_CTRL_STABLE_DT:
            self.dt = stable_timestep(self.T, self.const.GET_NU_FLOAT(self.RE), self.const.U, self.mesh.hmin())
        print "dt = {0:e}".format(self.dt)

        self.t = 0.0
        self.k = 0
        self.save_freq = int(1.0/(self.dt*self.const.LINEARIZED_CTRL_SAVE_PER_S))

        # read compress system
        names = ["M", "M_BOUNDARY_CTRL", "S", "R", "K", "G", "GT", "B", "C", "Kinf"]
        self.mat = {}
        for name in names:
            self.mat[name] = scio.mmread(const.ASSEMBLER_COMPRESS_CTRL_NAME_MTX(ref, name, RE))
        self.inner_nodes = np.loadtxt(const.ASSEMBLER_COMPRESS_CTRL_INNERNODES_DAT(ref, RE), dtype=np.int64)

        # log time, norm(u_delta), outputs, inputs
        self.logv = np.zeros((int(self.T/(self.dt*self.save_freq)+1),\
                              2+self.mat["C"].shape[0]+self.mat["B"].shape[1]))
        self.klog = 0
        self.kinfo = int(self.const.LINEARIZED_CTRL_INFO*(self.T/self.dt))

        # system sizes
        self.nv, self.np = self.mat["G"].shape
        self.ninner = self.inner_nodes.size

        # lift input matrix and feedback
        self.Kinfsys = np.vstack((self.mat["Kinf"], np.zeros((self.np, self.mat["Kinf"].shape[1]))))
        self.Bsys = -1*self.dt * np.vstack((self.mat["B"], np.zeros((self.np, self.mat["B"].shape[1]))))

        # visualisation
        self.u_dolfin = Function(self.V)
        self.u_dolfin.rename(*self.const.PVD_U_LABEL_NAME)
        self.u_file = File(const.LINEARIZED_CTRL_U_PVD(self.ref, self.RE))
        self.udelta_dolfin = Function(self.V)
        self.udelta_file = File(const.LINEARIZED_CTRL_U_DELTA_PVD(self.ref, self.RE))
        self.udelta_dolfin.rename(*self.const.PVD_U_LABEL_NAME)

        # define state
        self.uk_sys = np.take(self.u_stat.vector().array(), self.inner_nodes)
        self.uk_sys *= self.pertubationeps
        self.uk_uncps = np.zeros((self.V.dim(),))

        # build system matrices
        self.Asys = -self.mat["S"]-self.mat["R"]-self.mat["K"]-self.mat["M_BOUNDARY_CTRL"]
        u = scsp.hstack([self.mat["M"] - self.dt*self.Asys, self.dt*(-self.mat["G"])])
        l = scsp.hstack([self.dt * (-self.mat["GT"]), scsp.csr_matrix((self.np, self.np))])
        self.Msys_ode = scsp.vstack([u, l]).tocsc()
        self.Msys_lift = scsp.vstack([self.mat["M"], scsp.csr_matrix((self.np, self.nv))]).tocsr()

        # warning for higher refinements there are sometimes error change these parameters or ilu<->lu
        self.Msys_solver = scspli.splu(self.Msys_ode)

        # test function and form for assembleN
        self.w_test = TestFunction(self.V)
        self.N = inner(self.w_test, grad(self.u_dolfin) * self.u_dolfin)*dx
        self.N_uk_uk = np.zeros((self.np+self.ninner,))

        # build lu for smv
        temp = self.Msys_solver.solve(self.Bsys)
        temp = np.dot(self.Kinfsys.T, temp)
        temp = np.eye(temp.shape[0]) - temp
        lu, piv = scla.lu_factor(temp)
        self.smvlu_piv = (lu, piv)

    def assembleN(self, uk):
        # insert values at inner nodes
        self.uk_uncps[self.inner_nodes] = uk

        # update vector, assemble form  and compress vector
        self.u_dolfin.vector().set_local(self.uk_uncps)
        Nassemble = assemble(self.N)

        # compress to inner nodes and fill into upper block
        self.N_uk_uk[:self.ninner] = np.take(Nassemble.array(), self.inner_nodes)

    def uk_next(self):
        # compute new rhs (prediction)
        self.assembleN(self.uk_sys)
        Msys_lift_uk = self.Msys_lift * self.uk_sys
        rhs = Msys_lift_uk - self.dt*self.N_uk_uk

        if self.t < self.const.LINEARIZED_CTRL_START_CONTROLLING:
            # print "no feed"
            temp = self.Msys_solver.solve(rhs)
            # lift down vector
            self.uk_sys = temp[0:self.nv]
            # try to correct solution
            self.correction_no_ctrl(Msys_lift_uk)

        else:
            # perform sherman morrison woodbury
            temp = smw(self.smvlu_piv, self.Msys_solver, self.Bsys, self.Kinfsys, rhs)

            # lift down vector
            self.uk_sys = temp[0:self.nv]
            # try to correct solution
            self.correction_ctrl(Msys_lift_uk)

        self.k += 1
        self.t += self.dt

    def correction_ctrl(self, Msys_lift_uk):

        for i in range(self.const.LINEARIZED_CTRL_CORRECTION_STEPS):
            # assemble nonlinear right hand side term
            self.assembleN(self.uk_sys)
            rhs = -self.dt*self.N_uk_uk+Msys_lift_uk

            # perform sherman morrision woodbury
            ukpk_sys_correction = smw(self.smvlu_piv, self.Msys_solver, self.Bsys, self.Kinfsys, rhs)
            uk_sys_correction = ukpk_sys_correction[0:self.nv]


            if (i % self.const.LINEARIZED_CTRL_CORRECTION_RES_MOD) == 0:
                # compute residual for implicit euler scheme
                self.assembleN(uk_sys_correction)
                res = np.linalg.norm(self.Msys_ode.dot(ukpk_sys_correction)-\
                                     self.Bsys.dot(self.Kinfsys.T.dot(ukpk_sys_correction))+\
                                     self.dt*self.N_uk_uk-Msys_lift_uk)

                # show info
                if self.k % self.kinfo == 0:
                    print "Correction Step {0:d}\t||res (implicit euler)||={1:e}".format(i,  res)

                if np.isnan(res):
                    raise ValueError('nan during computation')

                if res < self.const.LINEARIZED_CTRL_CORRECTION_RES:
                    break

            self.uk_sys = uk_sys_correction

    def correction_no_ctrl(self, Msys_lift_uk):

        for i in range(self.const.LINEARIZED_CTRL_CORRECTION_STEPS):
            # assemble nonlinear right hand side term
            self.assembleN(self.uk_sys)
            ukpk_sys_correction = self.Msys_solver(-self.dt*self.N_uk_uk+Msys_lift_uk)
            uk_sys_correction = ukpk_sys_correction[0:self.nv]

            # compute residual for implicit euler scheme
            self.assembleN(uk_sys_correction)
            res = np.linalg.norm(self.Msys_ode*ukpk_sys_correction+self.dt*self.N_uk_uk-Msys_lift_uk)

            # show info
            if self.k % self.kinfo == 0:
                diff = np.linalg.norm(self.uk_sys-uk_sys_correction)
                print "Correction Step {0:d}\t||res (implicit euler)||={1:e}".format(i, res)

            self.uk_sys = uk_sys_correction
            if res < self.const.LINEARIZED_CTRL_CORRECTION_RES:
                break

    def log(self):
        # log time and u delta
        if self.k % self.save_freq == 0:
            self.uk_uncps[self.inner_nodes] = self.uk_sys
            self.u_dolfin.vector().set_local(self.uk_uncps)
            nrm = norm(self.u_dolfin, "L2", mesh=self.mesh)
            self.logv[self.klog, 0] = self.t
            self.logv[self.klog, 1] = nrm

            # log control
            uc = self.mat["Kinf"].T.dot(self.uk_sys)
            uc.ravel()
            self.logv[self.klog, 2:2+uc.size] = uc

            # log output
            uout = self.mat["C"].dot(self.uk_sys)
            uout.ravel()
            self.logv[self.klog, 2+uc.size:2+uc.size+uout.size] = uout
            self.klog += 1

            if np.isnan(nrm):
                print "nan during computation"
                return False

            if nrm < self.const.LINEARIZED_CTRL_RES:
                print "convergenced reached"
                return False

        return True

    def save(self):

        if self.k % self.save_freq == 0:
            # uncompress and save pvd for udelta
            self.uk_uncps[self.inner_nodes] = self.uk_sys
            self.udelta_dolfin.vector().set_local(self.uk_uncps)
            self.udelta_file << (self.udelta_dolfin, self.t)

            # add stationary and save pvd
            v = self.uk_uncps + self.u_stat.vector().array()
            self.u_dolfin.vector().set_local(v)
            self.u_file << (self.u_dolfin, self.t)

    def solve_ode(self):
        while self.t < (self.T + DOLFIN_EPS):
            self.save()

            if not self.log():
                break

            # print info
            if self.k % int(self.const.LINEARIZED_CTRL_INFO*(self.T/self.dt)) == 0:
                self.uk_uncps[self.inner_nodes] = self.uk_sys
                self.u_dolfin.vector().set_local(self.uk_uncps)
                nrm = norm(self.u_dolfin, "L2", mesh=self.mesh)
                print gettime(), "{0:.2f}%\t t={1:.3f}\t ||u_delta||={2:e}".format(self.t/self.T*100, self.t, nrm)

            self.uk_next()

    def save_log(self):
        # self.logv = self.logv[:self.k, :]
        np.savetxt(self.const.LINEARIZED_CTRL_LOG(self.ref, self.RE), self.logv)