# Script for running benchmark on cluster.
# Important: Must run combination of solver/problem to compile
# before submitting jobs on a higher refinement level.

# Submission is paused for 5 seconds to prevent small jobs from accessing the result.log file simultaneously.

# Kristian Valen-Sendstad, 2009

import os
from dolfin_utils.pjobs import submit
from time import *

#from pjobs import submit
if ".." not in os.environ['PYTHONPATH']:
    os.environ['PYTHONPATH'] = "..:" + os.environ['PYTHONPATH']
jobs = []

# Solvers and problems
solvers = ["chorin","css1","css2","ipcs","grpc"]
solvers = ["chorin"]
#solvers = ["css1"]
#solvers = ["css2"]
#solvers = ["ipcs"]
#solvers = ["grpc"]
#solvers = ["stat_newton"]


problems = ["drivencavity","channel","beltrami","cylinder","karman.karman"]
#problems = ["drivencavity"]
#problems = ["channel"]
#problems = ["beltrami"]
#problems = ["cylinder"]
#problems = ["karman.karman"]


# Number of refinement levels
refinements = [1]

# See output on screen
to_screen = True

# Max number of runtime days
days = 10

for k in range(len(refinements)):
        for j in range(len(problems)):
            jobs = []
            for i in range(len(solvers)):
                if to_screen==True:
                    #os.system("instant-clean")
                    print "python ns " + problems[j], solvers[i], "refinement_level=" +str(refinements[k]), ">&", problems[j] + "_" + solvers[i] + "_" + str(refinements[k])
                    jobs.append(("python ns " + problems[j] + " "+ solvers[i] + " " + "refinement_level=" +str(refinements[k]) + " "+ "debug=True" + ">&" + " " + problems[j]  + solvers[i] + "_" + str(refinements[k])))
                    os.system("python ns " + problems[j] + " "+ solvers[i] + " " + "refinement_level=" +str(refinements[k]) + " "+ "debug=True")


            #print ''
            #if problems[j] == "beltrami":
                #submit(jobs, nodes=1, ppn=8,  keep_environment=True, walltime=24*days, dryrun=False)
            #else:
                #submit(jobs, nodes=1, ppn=2,  keep_environment=True, walltime=24*days, dryrun=False)
            #print ''
            #print 'Submitted jobs:' , len(jobs)
            #print ''


            sleep(5)



# For larger jobs (More memory demanding jobs, might require all available memory on node =>
# args nodes=1, ppn=8 allocates entire node for one job.)
#submit(jobs, nodes=1 , ppn=8 , keep_environment=True, walltime=24*days, dryrun=False)
