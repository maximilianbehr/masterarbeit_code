from src.stationary.newton import Newton
import traceback


# set Reynoldsnumbers and refinements
REs = range(100, 5000, 100)
refs = [1, 2, 3, 4]



for ref in refs:
    REinitial = None
    for RE in REs:
        print "ref", ref
        print "RE", RE

        try:
            newton = Newton(ref, RE, REinitial)
            #newton = Newton(ref, RE, None)
            newton.solve()
            newton.save()
            REinitial = RE
        except Exception, e:
            traceback.print_exc()
            # reynolds number to large next refinement level
            break