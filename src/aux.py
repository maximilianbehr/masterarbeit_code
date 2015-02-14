import time
import os
from functools import wraps


def createdir(file):
    if not os.path.exists(os.path.dirname(file)):
        os.makedirs(os.path.dirname(file))

def gettime():
    return time.strftime("%d/%m/%Y %H:%M:%S")


PROF_DATA = {}
def profile(fn):
    @wraps(fn)
    def with_profiling(*args, **kwargs):
        start_time = time.time()

        ret = fn(*args, **kwargs)

        elapsed_time = time.time() - start_time
        fname = fn.__module__ +" "+fn.__name__
        if fname not in PROF_DATA:
            PROF_DATA[fname] = [0, []]
        PROF_DATA[fname][0] += 1
        PROF_DATA[fname][1].append(elapsed_time)

        return ret

    return with_profiling

def print_prof_data():
    for fname, data in PROF_DATA.items():
        max_time = max(data[1])
        avg_time = sum(data[1]) / len(data[1])
        print "Function %s called %d times. " % (fname, data[0]),
        print 'Execution time max: %.6f, average: %.6f' % (max_time, avg_time)

def clear_prof_data():
    global PROF_DATA
    PROF_DATA = {}