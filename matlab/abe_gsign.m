function [X,it] = abe_gsign(A,B,E)
%ABE_GSIGN
%
% Solve the anti-stable Algebraic Bernoulli Equation (ABE) 
%
% (*)    A' X E + E' X A - E' X B X E = 0
%   
% via the generalized Newton iteration for the sign function.
%
% INPUT: A, B, E the coefficient matrices in (*).
%        It is required that A, B, E are square, (A,E) must not 
%        have eigenvalues on the imaginary axis.        
%        NOTE: the assumption on the spectra is not checked!   
%
% OUTPUT: X, the solution of (*), and it, the number of iterations
%         used. 

% Peter Benner, 2008/03/16. 

narginchk(3,3)
if (size(A,1) ~= size(B,1)),  
  error('A and B must have same number of rows!'),  
end
if (size(A,2) ~= size(B,2)),  
  error('A and B must have same number of columns!'),  
end
n = size(A,1);
%
it      = 0;
maxit   = 50;
tol     = 10*n*sqrt(eps); 
%
Err         = 1;
onemore     = 0;
convergence = Err <= tol;

Enrm = norm(E,1);
[~,UE,~] = lu(E);
de = abs(diag(UE));
if any(de < n*eps*Enrm),  warning('E must not be singular.'),  end
de = de.^(1/n);
de = prod(de);   

while (it < maxit) && ((~convergence) || (convergence && (onemore < 2))),
  [L,U,P] = lu(A);
  Ainv    = E*(U\(L\P));
% Determinantel scaling
  da      = abs(diag(U));
  cs      = prod(da.^(1/(n)))/de; 
  A       = (A/cs + cs*Ainv*E)/2;
  B       = (B/cs + cs*Ainv*B*Ainv')/2;
  Err  = norm(A - E,1)/Enrm;
  it   = it + 1;
  fprintf('Step %i, conv. crit. = %d\n', it, Err)
  convergence = Err <= tol;
  if convergence,  onemore = onemore + 1;  end
end
%
X = ([B; E' - A']\[A/E + eye(n); zeros(n)]);
X = (X + X')/2;

if (it == maxit && Err > tol),
  fprintf('ABE_SIGN: No convergence in %d iterations.\n', maxit)
end
