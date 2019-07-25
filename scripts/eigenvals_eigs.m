%% EIGENVALS_EIGS Computing egenvalues with the use of the descriptor matrix and Matlab's Eigs function (ARPACK Arnoldi)
function eigenvals_eigs()
global dif_eqs sigma first_time SPsolver Transp V eigenvals W eigenvals2 Tot_RHScorr_tmr Tot_RHSsol_tmr Tot_Injsol_tmr Tot_decomp_tmr Tot_factor_schur_tmr NumberofSols Fac_tmr Sol_tmr numEig

%% Select the parameters for eigs() and the solver to be used

% sigma = 1i ;
n = numel(dif_eqs);
opts.isreal = 0 ; % sigma can be complex, thus A-sigmaB can be complex
opts.issym = 0 ;
opts.tol = 10e-12 ; % convergence tolerance
% numEig = 4; % number of eigenvalues to be computed
opts.p = min(3*numEig+1,n) ;

%  Sparse solver to be used. The options available are in eigs_solver.m
SPsolver = 'KLU'; % MatLU, KLU, PARDISO, DPS
% dd = zeros(numEig,1);
Tot_RHScorr_tmr=0; 
Tot_RHSsol_tmr=0; 
Tot_Injsol_tmr=0;
Tot_decomp_tmr=0;
Tot_factor_schur_tmr=0;
% NumberofSols=0;
Fac_tmr=0;
Sol_tmr=0;

%% Compute the Eigenvalues and the right eigenvectors

Transp = false;
first_time = 1;
right_tmr=tic;
[V_tmp,eigenvals_tmp] = eigs(@eigs_solver,n,numEig,sigma,opts);
right_tmr=toc(right_tmr);
V = [V V_tmp];
eigenvals = [eigenvals ; diag(eigenvals_tmp)];
% display(eigenvals);
% fprintf('Eigenvalues and right eigenvectors using eigs done in %.3f seconds. Took %d sparse solutions.\n',right_tmr,NumberofSols);

%% Compute the Eigenvalues and the left eigenvectors
%  This is done by using the transpose of the descriptor matrix in the
%  Arnoldi iterations
Transp = true;
first_time = 1;
NumberofSols = 0;
left_tmr=tic;
[W_tmp,eigenvals2_tmp] = eigs(@eigs_solver,n,numEig,sigma,opts);
left_tmr=toc(left_tmr);
eigenvals2 = [eigenvals2 ; diag(eigenvals2_tmp)];
W = [W W_tmp];
% fprintf('Eigenvalues and left eigenvectors using eigs done in %.3f seconds. Took %d sparse solutions.\n',left_tmr,NumberofSols);
% fprintf('Profiling:\n==========\n');
% fprintf('Solver initialization and factoriz. done in %.3f seconds.\n',Fac_tmr);
% fprintf('Solver solutions done in %.3f seconds.\n\n',Sol_tmr);
% if strcmp(SPsolver,'Decomposed')
%     fprintf('(P) Decomposition of system done in %.3f seconds.\n',Tot_decomp_tmr);
%     fprintf('Factorization of Schur-complement done in %.3f seconds.\n',Tot_factor_schur_tmr);
%     fprintf('(P) Correction of RHS for Schur-complement done in %.3f seconds.\n',Tot_RHScorr_tmr);
%     fprintf('Solution of Schur-complement done in %.3f seconds.\n',Tot_RHSsol_tmr);
%     fprintf('(P) Solution of injectors done in %.3f seconds.\n',Tot_Injsol_tmr);
%     fprintf('Parallel percentage is %.2f.\n',(Tot_decomp_tmr+Tot_RHScorr_tmr+Tot_Injsol_tmr)/(right_tmr+left_tmr)*100.0);
% end
%% Print the computed eigenvalues in descending real order
% [foo, sortidx] = sort(real(diag(eigenvals_tmp)),'descend');
% dd=diag(eigenvals_tmp);
% dd=dd(sortidx);
% fprintf('\nsigma= %4.2g + j%4.2g\nSorted eigenvalues (eigs):\n',real(sigma),imag(sigma));
% for i = 1:numEig
%     if imag(dd(i)) >= 0
%         fprintf('Eigenvalue %3d:\t%12.6f+%12.6fi\tf= %f\n',i,real(dd(i)),imag(dd(i)),abs(imag(dd(i)))/(2*pi));
%     end
% end

%% Plot the computed eigenvalues
% figure(fig)
% plot(real(dd),imag(dd),'bd')
% assignin('base', 'arnoldi_nordic', dd);
% assignin('base', 'arnoldi_expanded_nordic', dd);
% hold on % This is to put the eigenvalues computed by eig on the same diagram
end