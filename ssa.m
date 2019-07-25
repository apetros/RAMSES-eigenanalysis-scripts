%% SSA Receives the descriptor form of an eigenvalue problem and calculates the eigenvalues/eigenvectors with different methods
%  Usage:
%  ssa("jac_val.dat","jac_eqs.dat","jac_var.dat","jac_struc.dat",real_limit,damp_ratio)
%  ssa("jac_val.dat","jac_eqs.dat","jac_var.dat","jac_struc.dat",real_limit)
%  ssa("jac_val.dat","jac_eqs.dat","jac_var.dat","jac_struc.dat")
%
%  jac_val : has the values of the matrix in coordinate format
%  jac_eqs : has the description of the equations, mainly if they are
%  differential or algebraic
%  jac_val : has the description of the variables, mainly if they are
%  differential or algebraic
%  jac_struc : decomposed structure of power system
%  real_limit : the real number above which the eigenvalue is considered as
%  dominant (optional, default=-inf)
%  damp_ratio : Dumping ratio above which the eigenvalue is considered as
%  dominant (optional, default=1.0)
function ssa(jac_val,jac_eqs,jac_vars,jac_struc,real_limit,damp_ratio)

global sigma analysis numEig

%% Check and correct arguments
if nargin < 3
    error('You need to give the description of the matrix. Write "help ssa".');
elseif nargin > 6
    error('ssa requires at most 5 inputs. Write "help ssa".');
elseif nargin < 4
    jac_struc = '';
    real_limit = -inf;
    damp_ratio = 1.0;
elseif nargin < 5
    real_limit = -inf;
    damp_ratio = 1.0;
elseif nargin < 6
    damp_ratio = 1.0;
end

%% Initialize data structure
addpath('scripts') % Adds path to subfolder with scripts
init_tmr = tic;
init(jac_val, jac_eqs, jac_vars, jac_struc);
fprintf('Initialization done in %.3f seconds.\n\n', toc(init_tmr));

%% Calculate eigenvalues and eigenvectors using eigs (ARPACK Arnoldi method)
% Eigs_tmr=tic;
% analysis = 'IRA';
% numEig = 10 ;
% real_part = -0.1;
% while real_part > -1.0
%     freq = sqrt((real_part^2-damp_ratio^2*real_part^2)/damp_ratio^2)/(2.0*pi);
%     while freq < 2.0
%         sigma = real_part+2.0*pi*freq*1i;
% %         sigma = 0.0 + 6.28*1i;
% %         freq=6.28/(2.0*pi);
%         fprintf('Sigma= %.2f + i %.2f (freq= %.2f).\n',real(sigma),imag(sigma),freq);
%         eigenvals_eigs();
%         sigma = real_part-2.0*pi*freq*1i;
% %         sigma = -0.5 - 2.0*pi*0.5*1i;
%         fprintf('Sigma= %.2f + i %.2f (freq= %.2f).\n',real(sigma),imag(sigma),freq);
%         eigenvals_eigs();
%         freq = freq +0.1;
%     end
%     real_part = real_part - 0.1;
% end
% sigma = 0.1*1i;
% eigenvals_eigs();
% analyze_results(real_limit,damp_ratio);
% fprintf('Total time spent in eigs %.3f seconds.\n\n',toc(Eigs_tmr));

%% Calculate eigenvalues and eigenvectors using JDQR
% Eigs_tmr=tic;
% eigenvals_jdqr();
% fprintf('Total time spent in eigs %.3f seconds.\n\n',toc(Eigs_tmr));

%% Calculate eigenvalues and eigenvectors using eig (QZ method)
Eig_tmr=tic;
analysis = 'QZ';
eigenvals_eig();
fprintf('Total time spent in eig %.3f seconds.\n\n',toc(Eig_tmr));
analyze_results(real_limit, damp_ratio);
evalin('base', 'save(''modal_reduction'')');
evalin('base', 'savefig(''eigs'')');
loop_analysis();
% fprintf('\nExecute ssa to rerun the analysis.\n');
end
