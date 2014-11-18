from dolfin import *




def newton_solve(RE, refalg, reflevel):
    # Parameters
    nu = Constant(0.01)
    dt = 0.1
    idt = Constant(1. / dt)
    t_end = 10
    theta = 0.5  # Crank-Nicholson timestepping
    # theta   = 1   # Crank-Nicholson timestepping


    # Load problem_mesh and boundary function from file
    mesh = Mesh(KARMAN_REFN_MESH_FILE(refalg, reflevel))
    boundaryfunction = MeshFunction("size_t", mesh, KARMAN_REFN_BOUNDARY_FILE(refalg, reflevel))

    # Define function spaces (P2-P1)
    V = VectorFunctionSpace(mesh, "Lagrange", 2)
    Q = FunctionSpace(mesh, "Lagrange", 1)
    W = V * Q

    # define test functions
    (v, q) = TestFunctions(W)

    # define trial function
    w = Function(W)
    (u, p) = (as_vector((w[0], w[1])), w[2])

    # boundary parts Measure
    ds = Measure("ds")[boundaryfunction]

    # No-slip boundary condition for velocity
    noslip = Constant((0, 0))
    u_in = Expression(("(1-x[1])*x[1]*2", "0"))
    noslip_lower = DirichletBC(W.sub(0), noslip, GammaLower())
    noslip_upper = DirichletBC(W.sub(0), noslip, GammaUpper())
    noslip_ball = DirichletBC(W.sub(0), noslip, GammaBall())
    inflow = DirichletBC(W.sub(0), u_in, GammaLeft())
    bc = [noslip_upper, noslip_lower, noslip_ball, inflow]

    n = FacetNormal(mesh)
    I = Identity(V.cell().d)  # Identity tensor

    # current master step
    w = Function(W)
    (u, p) = (as_vector((w[0], w[1])), w[2])
    D = 0.5 * (grad(u) + grad(u).T)
    T = -p * I + 2.0 * nu * D

    # previous master step
    w0 = Function(W)
    (u0, p0) = (as_vector((w0[0], w0[1])), w0[2])
    D0 = 0.5 * (grad(u0) + grad(u0).T)
    T0 = -p * I + 2.0 * nu * D0

    # Define variational forms without master derivative in previous master
    F0_eq1 = (inner(T0, grad(v))) * dx + inner(grad(u0) * u0, v) * dx
    F0_eq2 = 0 * q * dx
    F0 = F0_eq1 + F0_eq2

    # variational form without master derivative in current master
    F1_eq1 = (inner(T, grad(v)) + inner(grad(u) * u, v)) * dx
    F1_eq2 = q * div(u) * dx
    F1 = F1_eq1 + F1_eq2

    # combine variational forms with master derivative
    #
    # dw/dt + F(t) = 0 is approximated as
    # (w-w0)/dt + (1-theta)*F(t0) + theta*F(t) = 0
    #
    F = idt * inner((u - u0), v) * dx + (1.0 - theta) * F0 + theta * F1

    # residual of strong Navier-Stokes
    r = idt * (u - u0) + theta * grad(u) * u + (1.0 - theta) * grad(u0) * u0 \
        - theta * div(T) - (1.0 - theta) * div(T0)


    # define Jacobian
    J = derivative(F, w)

    # Create files for storing solution
    u_pvd_file = File(INSTAT_NSE_NEWTON_VELOCITY_PVD_FILE(RE, refalg, reflevel))
    p_pvd_file = File(INSTAT_NSE_NEWTON_PRESSURE_PVD_FILE(RE, refalg, reflevel))

    # create variational problem and solver
    problem = NonlinearVariationalProblem(F, w, bc, J)
    solver = NonlinearVariationalSolver(problem)
    solver.parameters['newton_solver']['maximum_iterations'] = 20

    # Time-stepping
    t = dt
    while t < t_end:
        print "t =", t

        # Compute
        begin("Solving ....")
        solver.solve()
        end()

        # Extract solutions:
        (u, p) = w.split()

        # Plot
        #plot.py(u)

        # Save to file
        u_pvd_file << (u, t)
        p_pvd_file << (p, t)

        # Move to next master step
        w0.assign(w)
        t += dt


if __name__ == "__main__":
    newton_solve(1, "bisection", int(sys.argv[1]))

