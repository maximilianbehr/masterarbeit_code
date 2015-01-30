from src.lqr.bernoulli import Bernoulli
import traceback


refs = [3]
REs = range(400, 1000, 100)

for ref in refs:
    for RE in REs:
        try:
            print "ref={0:d} RE={1:d}".format(ref, RE)
            bernoulli = Bernoulli(ref, RE)
            bernoulli.solve()
            bernoulli.save()
        except Exception, e:
            traceback.print_exc()