clear; clc;
rng(2021);

sys    = LTISystem;
sys.Nx = 10; alpha = 0.8; rho = 2; actDens = 0.6; 
generate_dbl_stoch_chain(sys, rho, actDens, alpha);
sys.B1 = eye(sys.Nx);

% Params
params = MPCParams();
params.locality_ = 3;
params.tFIR_     = 2;

params.stateConsMtx_     = eye(sys.Nx);
params.stateUB_          = 4 * ones(sys.Nx, 1); % phase
params.stateUB_(2:2:end) = 10; % frequency
params.stateLB_          = -params.stateUB_;    
params.QSqrt_ = eye(sys.Nx);
params.RSqrt_ = eye(sys.Nu);

[Hset, hset] = terminal_set(sys, params);

% Since this matrix adds coupling, we need to keep track of who couples
% with who
% for i = 1:Nx
% rowstokeep{i} = find(H_terminal(:,i));
% H_ter_local{i} = H_terminal(rowstokeep{i},:); % keep only rows concerning subsystem i
% h_ter_local{i} = h_terminal(rowstokeep{i});
% neighbors{i} = find(~any(H_ter_local{i},1)==0); % keep track of neighbors of subsystem i
% H_ter_local{i} = H_ter_local{i}(:,neighbors{i}); % eliminate the rows with only zeros    
% end
