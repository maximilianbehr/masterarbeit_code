import time
import os
from functools import wraps
import gzip
import scipy.io as scio
import os

def write_matrix(name, mat, com):
    print "Write Matrix {0:s}".format(name)
    pid = os.getpid()
    scio.mmwrite("tmp.dat_{0:d}.mtx".format(pid), mat, comment=com)
    f_in = open("tmp.dat_{0:d}.mtx".format(pid), "rt")
    f_out = gzip.open(name, "wt", compresslevel=9)
    f_out.writelines(f_in)
    f_out.close()
    f_in.close()
    #os.remove("tmp.dat.mtx")

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
