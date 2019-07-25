%% Decomposes the power system into network VS injectors
function decomp_fact_system()

global Sp nbbus adf nbsync nbinj nbtwop bus_inj A Dt Bx By Ctx Cty LUDt LUA Tot_decomp_tmr Tot_factor_schur_tmr
decomp_tmr=tic;

A = cell(nbsync+nbinj+nbtwop,1);
Bx = cell(nbsync+nbinj+nbtwop,1);
By = cell(nbsync+nbinj+nbtwop,1);
Cx = cell(nbsync+nbinj+nbtwop,1);
Cy = cell(nbsync+nbinj+nbtwop,1);
Bt = cell(nbsync+nbinj+nbtwop,1);
Ctx = cell(nbsync+nbinj+nbtwop,1);
Cty = cell(nbsync+nbinj+nbtwop,1);
LUA = cell(nbsync+nbinj+nbtwop,1);
for j=1:nbsync+nbinj+nbtwop
    A{j}=Sp(2*nbbus+adf(j):2*nbbus+adf(j+1)-1,2*nbbus+adf(j):2*nbbus+adf(j+1)-1);
    LUA{j} = klu (A{j});
    
    Bx{j}=Sp(2*nbbus+adf(j):2*nbbus+adf(j+1)-1,2*bus_inj(j)-1);
    By{j}=Sp(2*nbbus+adf(j):2*nbbus+adf(j+1)-1,2*bus_inj(j));
    
    Cx{j}=Sp(2*bus_inj(j),2*nbbus+adf(j):2*nbbus+adf(j+1)-1);
    Cy{j}=Sp(2*bus_inj(j)-1,2*nbbus+adf(j):2*nbbus+adf(j+1)-1);
        
    Cty{j}=(A{j}'\Cy{j}')';
    Ctx{j}=(A{j}'\Cx{j}')';
    
    Bt{j}(1)=Cty{j}*Bx{j};
    Bt{j}(2)=Cty{j}*By{j};
    Bt{j}(3)=Ctx{j}*Bx{j};
    Bt{j}(4)=Ctx{j}*By{j};
end
N=double(2*nbbus);
Corr = sparse(N,N);
for j=1:nbsync+nbinj+nbtwop
    i=bus_inj(j);
    Corr(2*i-1,2*i-1) = Corr(2*i-1,2*i-1) + Bt{j}(1);
    Corr(2*i-1,2*i) = Corr(2*i-1,2*i) + Bt{j}(2);
    Corr(2*i,2*i-1) = Corr(2*i,2*i-1) + Bt{j}(3);
    Corr(2*i,2*i) = Corr(2*i,2*i) + Bt{j}(4);
end
Dt = Sp(1:2*nbbus,1:2*nbbus)-Corr;
Tot_decomp_tmr=Tot_decomp_tmr+toc(decomp_tmr);
factor_tmr=tic;
LUDt = klu (Dt);
Tot_factor_schur_tmr=Tot_factor_schur_tmr+toc(factor_tmr);
end
