from dolfin import *
import numpy as np
from src.aux import profile, gettime
from linearized_aux import stable_timestep


#set floating point error handling
np.seterr(all="raise", divide="raise", over="raise", under="raise", invalid="raise")

class LinearizedSimPETSC():

    def __init__(self, const, ref, RE):
        parameters["reorder_dofs_serial"] = False
        parameters["reorder_vertices_gps"] = False
        parameters["reorder_cells_gps"] = False

        # parameters
        self.ref = ref
        self.RE = RE
        self.nu = const.GET_NU(RE)
        self.const = const
        self.pertubationeps = self.const.LINEARIZED_SIM_PERTUBATIONEPS
        self.penalty_eps = const.ASSEMBLER_PENALTY_EPS
        self.dt = self.const.LINEARIZED_SIM_DT
        self.T = self.const.LINEARIZED_SIM_T

        # mesh and function spaces, stationary solution and dirichlet bcs
        self.mesh = Mesh(const.MESH_XML(ref))
        self.V = VectorFunctionSpace(self.mesh, const.V, const.V_DIM)
        self.Q = FunctionSpace(self.mesh, const.Q, const.Q_DIM)
        self.W = self.V*self.Q
        self.boundaryfunction = MeshFunction("size_t", self.mesh, const.BOUNDARY_XML(ref))
        self.u_stat = Function(self.V, const.STATIONARY_U_XML(ref, RE))
        self.w_stat = Function(self.W, const.STATIONARY_W_XML(ref, RE))

        # get zero dirichlet boundary conditions
        noslip = Constant((0.0, 0.0))
        self.bcups = [DirichletBC(self.W.sub(0), noslip, self.boundaryfunction, INDICES) \
                      for INDICES in self.const.ASSEMBLER_COMPRESS_SIM_DIRI_ZEROS]

        if const.LINEARIZED_SIM_STABLE_DT:
            # characterteristic velocity is chose as 1
            self.dt = stable_timestep(self.T, self.const.GET_NU_FLOAT(self.RE), 1, self.mesh.hmin())
        else:
            self.dt = self.const.LINEARIZED_SIM_DT
        print "dt = {0:e}".format(self.dt)

        self.t = 0.0
        self.k = 0
        self.save_freq = int(1.0/(self.dt*self.const.LINEARIZED_SIM_SAVE_PER_S))

        # log
        self.logv = np.zeros((int(self.T/(self.dt*self.save_freq)+1), 2))
        self.klog = 0
        self.kinfo = int(self.const.LINEARIZED_SIM_INFO*(self.T/self.dt))

        # build system matrices from weak formulation
        # trial and test functions
        (dudt, dpdt) = TrialFunctions(self.W)
        (u, p) = TrialFunctions(self.W)
        (w_test, p_test) = TestFunctions(self.W)
        ds = Measure("ds")[self.boundaryfunction]
        M = inner(dudt, w_test) * dx
        S = self.nu * inner(grad(u), grad(w_test)) * dx
        K = inner(grad(self.u_stat) * u, w_test) * dx
        R = inner(grad(u) * self.u_stat, w_test) * dx
        G = p * div(w_test) * dx
        GT = div(u) * p_test * dx

        A = -1*(S+R+K)
        self.Msys_ode = assemble(M-1*self.dt*A-1*self.dt*G-1*self.dt*GT)
        self.Msys_lift = assemble(M)

        # set initial value
        self.up_k = Function(self.W)
        self.up = project(self.const.LINEARIZED_SIM_PERTUBATIONEPS*self.w_stat, self.W)

        # apply dirichlet bcs dont pertubate dirichlet nodes
        [bcup.apply(self.up.vector()) for bcup in self.bcups]
        [bcup.apply(self.Msys_lift) for bcup in self.bcups]
        [bcup.apply(self.Msys_ode) for bcup in self.bcups]

        # for visualisation
        self.u_file = File(const.LINEARIZED_SIM_U_PVD(self.ref, self.RE))
        self.u = Function(self.V)
        self.p = Function(self.Q)
        self.udelta_file = File(const.LINEARIZED_SIM_U_DELTA_PVD(self.ref, self.RE))
        self.udelta = Function(self.V)
        self.pdelta = Function(self.Q)

        # setup lu solver
        self.solver = LUSolver()
        self.solver.set_operator(self.Msys_ode)
        self.solver.parameters["reuse_factorization"] = True

    @profile
    def assembleN(self):

        (w_test, p_test) = TestFunctions(self.W)
        (u, p) = self.up.split()
        self.N = assemble(inner(w_test, grad(u) * u)*dx)
        [bcup.apply(self.N) for bcup in self.bcups]

    @profile
    def uk_next(self):

        self.assembleN()
        Msys_liftup = self.Msys_lift*self.up.vector()
        rhs = self.Msys_lift*self.up.vector()-self.dt*self.N
        self.solver.solve(self.up_k.vector(), rhs)
        # print norm(self.Msys_ode*self.up_k.vector()-rhs, 'l2')
        self.up.assign(self.up_k)
        # try to correct solution
        self.correction(Msys_liftup)
        self.k += 1
        self.t += self.dt

    def correction(self, Msys_liftup):
        # correct predicted solution
        for i in range(self.const.LINEARIZED_SIM_CORRECTION_STEPS):
            # assemble nonlinear right hand side term
            self.assembleN()
            self.solver.solve(self.up_k.vector(), -self.dt*self.N+Msys_liftup)
            self.up.assign(self.up_k)

            # compute residual for implicit euler scheme only in every third step to improve perfomance
            if (i % self.const.LINEARIZED_SIM_CORRECTION_RES_MOD) == 0:
                self.assembleN()
                residual= self.Msys_ode*self.up.vector()+self.dt*self.N-Msys_liftup
                res = np.linalg.norm(residual.array())

                # show info
                if self.k % self.kinfo == 0:
                    print "Correction Step {0:d}\t||res (implicit euler)||={1:e}".format(i,  res)

                if np.isnan(res):
                    raise ValueError('nan during computation')

                if res < self.const.LINEARIZED_SIM_CORRECTION_RES:
                    break

    def log(self):
        if self.k % self.save_freq == 0:
            # nrm = np.linalg.norm(self.uk_sys)
            (self.u, self.p) = self.up.split()
            nrm = norm(self.u, "L2", mesh=self.mesh)
            self.logv[self.klog, 0] = self.t
            self.logv[self.klog, 1] = nrm

            self.klog += 1

            if np.isnan(nrm):
                print "nan during computation"
                return False

            if nrm < self.const.LINEARIZED_SIM_RES:
                print "convergenced reached"
                return False

        return True


    def save(self):
        if self.k % self.save_freq == 0:

            # uncompress and save pvd for udelta
            (self.udelta, self.pdelta) = self.up.split(deepcopy=True)
            self.udelta_file << (self.udelta, self.t)

            # add stationary and save pvd
            v = self.udelta.vector().array() + self.u_stat.vector().array()
            self.u.vector().set_local(v)
            self.u_file << (self.u, self.t)

    def solve_ode(self):
        while self.t < (self.T + DOLFIN_EPS):
            self.save()

            if not self.log():
                break

            # print info
            if self.k % self.kinfo == 0:
                (self.u, self.p) = self.up.split()
                print gettime(), "{0:.2f}%\t t={1:.3f}\t||u_delta||={2:e}".format(self.t/self.T*100, self.t, norm(self.u,'L2', mesh=self.mesh))

            self.uk_next()

    def save_log(self):
        # self.logv = self.logv[:self.k, :]
        np.savetxt(self.const.LINEARIZED_SIM_LOG(self.ref, self.RE), self.logv)



















