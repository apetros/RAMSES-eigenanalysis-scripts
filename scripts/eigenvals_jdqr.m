%% EIGENVALS_JDQR Computing egenvalues with the use of the state matrix and JD function
function eigenvals_jdqr()
global dif_eqs alg_eqs S Jdyn_sp dif_states alg_states V eigenvals W Jdyn Jdyn_sp

calc_Jdyn();

%% Compute the eigenvalues, right (V), and left (W) eigenvectors
Eig_calc_tmr=tic;
method='jdqr';
nselect=5;
SHOW=1;
sigma=-1.0+6i;
tol=1.e-1;
TS=0;
MaxIt=1000;
options=struct('Tol',tol,'Disp',SHOW,'TestSpace',TS,'MaxIt',MaxIt); 
[Xeig,Lambda]=feval(method,Jdyn_sp,nselect,sigma,options); 
nrm=norm(Jdyn_sp*Xeig-Xeig*Lambda);
if nrm> tol
   fprintf('\nnorm(%s*X-X*Lambda): %0.5g\n','Jdyn',nrm), pause
end
fprintf('Computation of eigenvalues, left, and right eigenvectors with JDQR done in %.3f seconds.\n',toc(Eig_calc_tmr));
end
