function [Feed0, Feed1]=initial_feedback_karman(Mmtx, Smtx,Rmtx,Kmtx,Mlowermtx, Muppermtx, Gmtx, Bmtx,Cmtx)
% Compute the initial feedback with projection method 
% H. Weichelt, J. Saak, P. Kürschner September 2011
%
% Feed0 -> (A-B*Feed0') is stable
% Feed1 -> (A-Feed1*C)  is stable

%%build full matrices

M = mmread(Mmtx);
S = mmread(Smtx);
R = mmread(Rmtx);
K = mmread(Kmtx);
Mlower = mmread(Mlowermtx);
Mupper = mmread(Muppermtx);
G = mmread(Gmtx);
B = mmread(Bmtx);
C = mmread(Cmtx);


[nv,np] = size(G);
fullA = [-S-R-K-Mlower-Mupper,G;G',zeros(np,np)];
%fullA = [-S-R-K,G;G',zeros(np,np)];
fullM = [M,zeros(nv,np);zeros(np,nv),zeros(np,np)];
%fullM = [M,-0.02*G;-0.02*G',zeros(np,np)];

% compute right and left eigenvectors
%neigs = ceil(size(fullA,1)/3);
%fprintf('compute %d smallest magnitude eigenvalues\n',neigs);
tic;
[vr,dr]=eigs(fullA, fullM, 250,'SM');
%[vr,dr]=eigs(fullA, fullM, neigs,'SM');
Ir=find(real(diag(dr))>0);
toc

tic;
[vl,dl]=eigs(fullA',fullM',250,'SM');
%[vl,dl]=eigs(fullA',fullM',neigs,'SM');
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

end
