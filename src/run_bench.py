import os

if ".." not in os.environ['PYTHONPATH']:
    os.environ['PYTHONPATH'] = "..:" + os.environ['PYTHONPATH']



#test problems from nsbench
instant_clean = False
refinements = [3]
problems = ["beltrami", "drivencavity"]
solvers = ["chorin", "css1", "css2", "ipcs"]
for refinement in refinements:
    for problem in problems:
        for solver in solvers:
            if instant_clean:
                os.system("instant-clean")
            print "python ns.py {0} {1} refinement_level={2} debug=True".format(problem, solver, str(refinement))
            os.system("python ns.py {0} {1} refinement_level={2} debug=True".format(problem, solver, str(refinement)))

#karman
instant_clean = False
REs = [100, 200, 300, 400, 500]
refinements = [2, 3, 4, 5]
problems = ["karman"]
solvers = ["chorin", "css1", "css2", "ipcs", "stat_newton"]
for RE in REs:
    for refinement in refinements:
        for problem in problems:
            for solver in solvers:
                if instant_clean:
                    os.system("instant-clean")
                print "python ns.py {0} {1} refinement_level={2} RE={3} debug=True".format(problem, solver, str(refinement),RE)
                os.system("python ns.py {0} {1} refinement_level={2} RE={3} debug=True".format(problem, solver, str(refinement),RE))
