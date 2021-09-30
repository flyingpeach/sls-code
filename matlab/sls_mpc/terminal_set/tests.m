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

[HTerminal, hTerminal] = terminal_set(sys, params);
sanity_check_terminal_coupling(sys, params.locality_+2, HTerminal);

params.terminal_H_ = HTerminal;
params.terminal_h_ = hTerminal;
