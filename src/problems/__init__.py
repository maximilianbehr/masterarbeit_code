# List of problems
problems = ["drivencavity", "beltrami", "karman"]

# Wrapper problem classes
def Problem(name, options):
    """Return problem instance for given problem name"""
    print name
    exec ("from %s import Problem as NamedProblem" % name)
    return NamedProblem(options)
