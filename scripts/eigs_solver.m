%% EIGS_SOLVER Implements y = A\x, needed for the Arnoldi method (eigs) using different solvers
function y = eigs_solver(x)
global dif_eqs alg_eqs myL myU sigma first_time S SLocal Sp LU SPsolver verbose info Transp NumberofSols Fac_tmr Sol_tmr

b = sparse(dif_eqs,ones(numel(dif_eqs),1),x,numel(dif_eqs)+numel(alg_eqs),1);
y = zeros(numel(dif_eqs),1);
NumberofSols=NumberofSols+1;

if first_time == 1
    eigs_solver_init=tic;
    if Transp
        SLocal = S';
    else
        SLocal = S;
    end
    
    if isnumeric(sigma) 
        Es = sparse(dif_eqs,dif_eqs,sigma*ones(1,numel(dif_eqs)),numel(dif_eqs)+numel(alg_eqs),numel(dif_eqs)+numel(alg_eqs)); % sigma*E
        Sp=SLocal-Es;
    else % If sigma is not a complex number, it means it's LM, SM, LR, etc.
        Sp=SLocal;
    end
    
    if strcmp(SPsolver,'MatLU') % Matlab LU and factorize
        [myL,myU]=lu(Sp); 
    elseif strcmp(SPsolver,'KLU')% KLU LU and factorize
        LU = klu (Sp); 
    elseif strcmp(SPsolver,'PARDISO')
        verbose = false ;
        info = pardisoinit(13,0);
        info.iparm(3) = 1; % number of threads to be used by pardiso
        info = pardisoreorder(Sp,info,verbose);
%         fprintf('The factors have %d nonzero entries.\n',info.iparm(18));
        info = pardisofactor(Sp,info,verbose); % Compute the numeric factorization.
    elseif strcmp(SPsolver,'DPS') % Decomposed parallel solver
        decomp_fact_system();
    else
        fprintf('Initialization of %s solver done in %.3f seconds.\n',SPsolver,toc(eigs_solver_init));
    end
    
    Fac_tmr=Fac_tmr+toc(eigs_solver_init);
    first_time = 0;
end

eigs_solver_sol=tic;
if strcmp(SPsolver,'MatLU')
    yy=myU\(myL\b); % Solve the descriptor matrix
elseif strcmp(SPsolver,'KLU')
    b = full(b); % KLU doesn't accept sparse RHS
    yy = klu(LU,'\',b) ;
elseif strcmp(SPsolver,'Mat')
    yy = Sp \ b;
elseif strcmp(SPsolver,'PARDISO')
    b = full(b); % PARDISO doesn't accept sparse RHS
    [yy info] = pardisosolve(Sp,b,info,verbose);
    % fprintf('PARDISO performed %d iterative refinement steps.\n',info.iparm(7));
elseif strcmp(SPsolver,'DPS')
    yy = solve_decomposed(b);
end

for idx = 1:numel(dif_eqs) % Copy back the states
    y(idx,1)=yy(dif_eqs(idx),1);
end

Sol_tmr=Sol_tmr+toc(eigs_solver_sol);
end