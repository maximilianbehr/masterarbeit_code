import numpy
import scipy

raise NotImplementedError('The ABE GSIGN Function iteration is not fully implemented')

A = -8*numpy.matrix(numpy.eye(6,6))
B = 3*numpy.matrix(numpy.eye(6,6))
E = 0.5*numpy.matrix(numpy.eye(6,6))


if not A.shape==B.shape:
    raise ValueError('A and B must have the same number of rows and cols!')

n = A.shape[0]
#parameters for iteration
it = 0
maxit = 50
eps = numpy.finfo(numpy.float).eps
tol = 10*n*numpy.sqrt(eps)
Err = 1
onemore = 0
convergence = Err <= tol

Enrm = numpy.linalg.norm(E, 1)
PE, LE, UE = scipy.linalg.lu(E)
de = numpy.abs(numpy.diag(UE))
if numpy.any(de < n*eps*Enrm):
    raise RuntimeWarning('E must not be singular.')
de = de**(1.0/n)
de = numpy.prod(de)


while (it < maxit) and ((not convergence) or (convergence and (onemore < 2))):
    P, L, U = scipy.linalg.lu(A)
    Ainv = numpy.dot(E, numpy.linalg.solve(U, numpy.linalg.solve(L, P)))
    #determinant scaling
    da = numpy.abs(numpy.diag(U))
    cs = numpy.float(numpy.prod(da**(1.0/n)))/de
    A = (A/cs + cs*Ainv*E)/2
    B = (B/cs + cs*Ainv*B*Ainv.T)/2;
    Err = numpy.linalg.norm(A-E, 1)/Enrm
    it += 1
    print "Step {0:d}, conv. crit. = {1:e}".format(it, Err)
    if convergence:
        onemore += 1
