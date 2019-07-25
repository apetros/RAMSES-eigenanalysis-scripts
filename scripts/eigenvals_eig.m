%% EIGENVALS_EIG Computing egenvalues with the use of the state matrix and Matlab's eig function
function eigenvals_eig()
global V eigenvals W Jdyn eigenvals2

V = [];
eigenvals = [];
eigenvals2 = [];
W = [];

calc_Jdyn();

%% Compute the eigenvalues, right (V), and left (W) eigenvectors
Eig_calc_tmr=tic;
[V, eigenvals, W] = eig(Jdyn);
eigenvals = diag(eigenvals);
fprintf('Computation of eigenvalues, left, and right eigenvectors with QZ done in %.3f seconds.\n',toc(Eig_calc_tmr));
end
