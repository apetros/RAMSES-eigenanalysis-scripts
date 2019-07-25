function calc_Jdyn()
global dif_eqs alg_eqs S dif_states alg_states Jdyn Jdyn_sp

Jdyn_calc_tmr=tic;
%% Create dynamic Jacobian
fx = S(dif_eqs,dif_states);
fy = S(dif_eqs,alg_states);
gx = S(alg_eqs,dif_states);
gy = S(alg_eqs,alg_states);

gygx=sparse(numel(alg_states),numel(dif_states));

for i = 1:numel(dif_states)
    gygx(:,i)=gy\gx(:,i);
end

Jdyn_sp = fx - fy * gygx ;
Jdyn = full(Jdyn_sp);
assignin('base', 'Jdyn_sp', Jdyn_sp);
assignin('base', 'Jdyn', Jdyn);
fprintf('Calculation of Jdyn done in %.3f seconds.\n',toc(Jdyn_calc_tmr));

