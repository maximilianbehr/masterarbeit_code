import src.bws.bws_const as const

import traceback


# set Reynoldsnumbers and refinements
REs = range(100, 1000, 100)
refs = [3]

for shear in [0, 1, 2, 3]:
#for shear in [3]:
    const.OUTPUTDIR_NAME = "results_bws_{0:d}".format(shear)
    from src.bws.mesh.bws import MeshBuilder
    from src.bws.newton_bws import Newton
    meshbuilder = MeshBuilder(const)
    for ref in range(max(refs)):
        meshbuilder.refine()
        meshbuilder.save()

    for ref in refs:
        REinitial = None
        for RE in REs:
            print "ref={0:d}\t RE={1:d}".format(ref,RE)
            try:
                newton = Newton(const, ref, RE, REinitial, shear)
                newton.solve()
                newton.save()
                REinitial = RE
            except Exception, e:
                traceback.print_exc()
                # reynolds number to large next refinement level
                break