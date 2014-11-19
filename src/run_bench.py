
import os

if ".." not in os.environ['PYTHONPATH']:
    os.environ['PYTHONPATH'] = "..:" + os.environ['PYTHONPATH']

# Solvers and problems
#solvers = ["chorin","css1","css2","ipcs", "stat_newton"]
#solvers = ["chorin","css1","css2","ipcs"]
#solvers = ["chorin"]
#solvers = ["css1"]
#solvers = ["css2"]
solvers = ["ipcs"]
#solvers = ["stat_newton"]


#problems = ["drivencavity"]
#problems = ["drivencavity"]
#problems = ["beltrami", "drivencavity"]
problems = ["karman"]


# Number of refinement levels
refinements = [3]

#clear instant cache
instant_clean = False

for refinement in refinements:
    for problem in problems:
        for solver in solvers:
            if instant_clean:
                os.system("instant-clean")
            print "python ns.py {0} {1} refinement_level={2} debug=True".format(problem, solver, str(refinement))

            os.system("python ns.py {0} {1} refinement_level={2} debug=True".format(problem, solver, str(refinement)))
