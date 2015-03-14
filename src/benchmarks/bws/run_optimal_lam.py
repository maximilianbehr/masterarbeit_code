# -*- coding: utf-8 -*-
import src.benchmarks.bws.bws_const as const
from src.benchmarks.aux import *
import numpy as np
from dolfin import *


const.OUTPUTDIR_NAME = "bws_optimal_lam"

def STATIONARY_U_PVD(ref, RE, lam, const):
    return os.path.join(const.OUTPUTDIR(), const.STATIONARY_DIR,\
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "u_{0:3f}.pvd".format(lam))

def STATIONARY_U_XML(ref, RE, lam, const):
    return os.path.join(const.OUTPUTDIR(), const.STATIONARY_DIR,\
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "u_{0:3f}.xml.gz".format(lam))

def STATIONARY_W_XML(ref, RE, lam, const):
    return os.path.join(const.OUTPUTDIR(), const.STATIONARY_DIR,\
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "w_{0:3f}.xml.gz".format(lam))

def STATIONARY_RESULTS(ref, RE, const):
    return os.path.join(const.OUTPUTDIR(), const.STATIONARY_DIR,\
                        "ref_{0:d}".format(ref), "RE_{0:d}".format(RE), "results.txt")


if __name__ == "__main__":
    # set Reynoldsnumbers and refinements and Parameters
    REs = range(50, 1500, 50)
    refs = [1]
    dtlam = 0.1
    lamstart = -0.75
    lamend = 0.75
    lams = np.arange(lamstart, lamend, dtlam)
    lams = np.append(lams, lamend)

    # build all meshes boundary functions ect.
    print "----------build mesh----------------"
    build_mesh(const, refs)
    print "----------finished build mesh-------"

    for ref in refs:
        # load mesh and shear function
        mesh = Mesh(const.MESH_XML(refs[0]))
        shear = MeshFunction("size_t", mesh, const.SHEAR_XML(refs[0], 3))
        ds = Measure("ds")[shear]
        V = VectorFunctionSpace(mesh, const.V, const.V_DIM)

        REold = None
        for RE in REs:
            results = np.zeros((2, lams.shape[0]))
            idx = 0

            for lam in lams:
                # compute stationary solution for given lambda
                const.STATIONARY_LAM = lam

                # load initial w if available and solve
                print "ref={0:d} RE={1:d} lam={2:f}".format(ref, RE, lam)
                try:
                    if REold:
                        initialw = STATIONARY_W_XML(ref, REold, lam, const)
                        if os.path.isfile(initialw):
                            newton = Newton(const, ref, RE, initialw)
                        else:
                            newton = Newton(const, ref, RE)
                    else:
                        newton = Newton(const, ref, RE)
                    newton.solve()
                except:
                    continue

                # compute cost
                w = newton.w
                u = interpolate(newton.u, V)
                dux_dy = grad(u[0])[1]
                fcost = assemble(dux_dy*(dux_dy-abs(dux_dy))*ds(const.GAMMASHEAR3_INDICES))

                # store result in array and save xml and pvd
                results[0, idx] = lam
                results[1, idx] = fcost
                File(STATIONARY_U_PVD(ref, RE, lam, const)) << u
                File(STATIONARY_U_XML(ref, RE, lam, const)) << u
                File(STATIONARY_W_XML(ref, RE, lam, const)) << w

                # increment idx
                idx += 1

                # print results
                print "------------LAMBDA={0:f} COST={1:e}------------".format(lam, fcost)

            # save results array
            np.savetxt(STATIONARY_RESULTS(ref, RE, const), results.T, fmt="%10.5f\t%e")

            # print minimum
            print "Minimum for RE={0:d} and ref={1:d}".format(RE, ref)

            print "Best 5 Values"
            results = results[:, results[1, :] > 0]
            print results[:, results[1, :].argsort()][:, 0:5]

            # set RE initial
            REold = RE

