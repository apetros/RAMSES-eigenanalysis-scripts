%% INIT Reads the descriptor matrix in raw format and creates the internal data structure
function init(jac_val, jac_eqs, jac_vars, jac_struc)
global dif_eqs alg_eqs dif_states alg_states raw_vars fig S nbbus adf nbsync nbinj nbtwop bus_inj V eigenvals W eigenvals2

%% Read values file
raw_val = importdata(jac_val) ;
S = spconvert(raw_val) ;

%% Read structure file
if ~isempty(jac_struc)
    fid = fopen (jac_struc);
    raw_struc = textscan(fid, '%d %d %d');
    fclose(fid);
end

nbbus  = raw_struc{1}(1);
nbsync = raw_struc{1}(2);
nbinj  = raw_struc{2}(2);
nbtwop = raw_struc{3}(2);

adf = zeros(nbsync+nbinj+nbtwop+1,1);
bus_inj = zeros(nbsync+nbinj+nbtwop+1,1);
for idx = 3:numel(raw_struc{2})
    adf(idx-2) = raw_struc{2}(idx);
    bus_inj(idx-2) = raw_struc{3}(idx);
end

%% Read equations file
fid = fopen (jac_eqs);
raw_eqs = textscan(fid, '%d %s %s %s %s %d');
fclose(fid);

%% Read variables file
fid = fopen (jac_vars);
raw_vars = textscan(fid, '%d %s %s %s %s');
fclose(fid);

dif_states = []; % saves the places of differential states
alg_states = []; % saves the places of algebraic states

for idx = 1:numel(raw_vars{2})
    if char(raw_vars{2}(idx)) == 'd'
%         disp(strcat('element: ',num2str(idx),' is differential'))
        dif_states = [dif_states idx];
    else
        alg_states = [alg_states idx];
    end
end

dif_eqs=[]; % saves the places of differential equations
alg_eqs=[]; % saves the places of algebraic equations

gamma = sparse(numel(raw_eqs{6}),1);
for idx = 1:numel(raw_eqs{6})
    gamma(idx)=raw_eqs{6}(idx);
end

assignin('base', 'gamma', gamma);

idx = 1;
while idx <= numel(raw_eqs{6})
    if (gamma(idx) ~= 0) && (gamma(idx) ~= idx)
        idx2 = gamma(idx);
%         disp(strcat('swapping: ',num2str(idx),' with ', num2str(idx2)))
        S([idx,raw_eqs{6}(idx)],:) = S([raw_eqs{6}(idx),idx],:);
        gamma([idx,idx2]) = gamma([idx2,idx]);
        idx = idx - 1;
    end
    idx = idx +1;
end

assignin('base', 'S', S);
    
for idx = 1:numel(raw_eqs{6})
    if gamma(idx) > 0
%         disp(strcat('element: ',num2str(idx),' is differential'))
        dif_eqs = [dif_eqs idx];
    else
        alg_eqs = [alg_eqs idx];
    end
end

V = [];
eigenvals = [];
eigenvals2 = [];
W = [];

Cdelta = sparse(length(dif_states),1);
Comega = sparse(length(dif_states),1);
B = sparse(length(dif_states),1);
for k = 1:length(dif_states)
    if strcmp(char(raw_vars{5}(dif_states(k))),'omega') 
        Cdelta(k)=1;
    elseif strcmp(char(raw_vars{5}(dif_states(k))),'delta')
        Comega(k)=1;
    elseif strcmp(char(raw_vars{5}(dif_states(k))),'vf')
        B(k)=1;
    end
end
assignin('base', 'Cdelta', Cdelta);
assignin('base', 'Comega', Comega);
assignin('base', 'B', B);

close all
hold on
fig = 1;

fprintf('Number of diff-alg states=%d\n', numel(dif_states) + numel(alg_states));
fprintf('Number of diff states=%d\n', numel(dif_states));
end