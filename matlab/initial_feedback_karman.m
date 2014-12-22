function [Feed0, Feed1]=initial_feedback_cylinderwake(Mmtx, Amtx, Gmtx, Bmtx,Cmtx)
% Compute the initial feedback with projection method 
% H. Weichelt, J. Saak, P. Kürschner September 2011
%
% Feed0 -> (A-B*Feed0') is stable
% Feed1 -> (A-Feed1*C)  is stable

%%build full matrices

M = mmread(Mmtx);
A = mmread(Amtx);
G = mmread(Gmtx);
B = mmread(Bmtx);
C = mmread(Cmtx);


[nv,np] = size(J);
fullA = [A,G;G',zeros(np,np)];
fullM = [M,zeros(nv,np);zeros(np,nv),zeros(np,np)];


% compute right and left eigenvectors
tic;
[vr,dr]=eigs(fullA, fullM, 1000,'SM');
Ir=find(real(diag(dr))>0);
toc

tic;
[vl,dl]=eigs(fullA',fullM',1000,'SM');
Il=find(real(diag(dl))>0);
toc

fprintf(1,'Found %d,%d instable eigenvalues.\n',size(Il,1),size(Ir,1));

if(size(Il,1)==0 && size(Ir,1)==0)
    fprintf(1,'No instable eigenvalues detected, no need for initial feedback!\n');
    Feed0=[];
    Feed1=[];
    return;
end

% Sort instabel EW and EV
LEV=vl(:,Il);
REV=vr(:,Ir);

% Projection
tA = LEV'*(fullA * REV);
tE = LEV'*(fullM * REV); % Identitaet nach Konstruktion
tB = LEV'* [B;sparse(np,size(B,2))];
tC = [C sparse(size(C,1),np)] * REV;

% Solve algebraic Bernoulli equation on projeted space
[XK,it] = abe_gsign(tA,tB*tB',tE);
fprintf(1,'ABE solved within %d steps.\n',it);

% Build initial feedback for nontransposed A
K0 = (([B;sparse(np,size(B,2))]'* LEV) * XK) * (LEV' * fullM);

% Solve algebraic Bernoulli equation on projeted space
[XL,it] = abe_gsign(tA',tC'*tC,tE');
fprintf(1,'ABE solved within %d steps.\n',it);

% Build initial feedback for transposed A
L0 = (([C sparse(size(C,1),np)] * REV) * XL) * (REV' * fullM');

Feed0=real(K0(:,1:nv))';
Feed1=real(L0(:,1:nv))';

mmwrite('Feed0.mtx',Feed0);
mmwrite('Feed1.mtx',Feed1);

end
