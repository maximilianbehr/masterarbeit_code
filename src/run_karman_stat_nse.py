from src.stationary.newton import Newton


# set Reynoldsnumbers and refinements
REs = range(100, 3000, 100)
refs = [4]

for ref in refs:
    REinitial = None
    for RE in REs:
        print "ref", ref
        print "RE", RE

        newton = Newton(ref, RE, REinitial)
        newton.solve()
        newton.save()

        REinitial = RE