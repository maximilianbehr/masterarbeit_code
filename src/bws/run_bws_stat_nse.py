import src.bws.bws_const as const
from src.stationary.newton import Newton
import traceback


# set Reynoldsnumbers and refinements
REs = range(100, 1100, 100)
refs = [3]


for ref in refs:
    REinitial = None
    for RE in REs:
        print "ref", ref
        print "RE", RE

        try:
            newton = Newton(const, ref, RE, REinitial)
            newton.solve()
            newton.save()
            REinitial = RE
        except Exception, e:
            traceback.print_exc()
            # reynolds number to large next refinement level
            break