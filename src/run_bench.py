import os

if ".." not in os.environ['PYTHONPATH']:
    os.environ['PYTHONPATH'] = "..:" + os.environ['PYTHONPATH']

# Solvers and problems
solvers = ["chorin", "css1", "css2", "ipcs", "stat_newton"]
#solvers = ["chorin","css1","css2","ipcs"]
#solvers = ["chorin"]
#solvers = ["css1"]
#solvers = ["css2"]
#solvers = ["ipcs"]
#solvers = ["stat_newton"]


#problems = ["drivencavity"]
problems = ["beltrami", "drivencavity", "karman"]
#problems = ["karman"]


# Number of refinement levels
refinements = [2,3,4,5]
#refinements = [2]


#Reynoldszahl
REs = [100,200,300,400,500]
#REs = [100]


#clear instant cache
instant_clean = False

for RE in REs:
    for refinement in refinements:
        for problem in problems:
            for solver in solvers:
                if instant_clean:
                    os.system("instant-clean")

                #use only stat_newton for karman
                if not solver == "stat_newton" or problem=="karman":
                    print "python ns.py {0} {1} refinement_level={2} RE={3} debug=True".format(problem, solver, str(refinement),RE)
                    os.system("python ns.py {0} {1} refinement_level={2} RE={3} debug=True".format(problem, solver, str(refinement),RE))

