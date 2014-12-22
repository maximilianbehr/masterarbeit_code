function ploteigs(eigmtx,eigeps)

close all, clc, clearvars -except eigmtx eigeps
eigenvalues = mmread(eigmtx);

stable_eigenvalues = eigenvalues(real(eigenvalues)<0);
unstable_eigenvalues = eigenvalues(real(eigenvalues)>=0);


semilogx(real(stable_eigenvalues),imag(stable_eigenvalues),'rx')
hold on
%semilogx(real(unstable_eigenvalues),imag(unstable_eigenvalues),'bx')
semilogx(-real(stable_eigenvalues),2*imag(stable_eigenvalues),'bx')
saveas(gcf,eigeps,'epsc')
%close all, clc, clear all






