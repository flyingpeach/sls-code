clear all; close all; clc;

%% Setup plant + parameters
rng(421);
sys    = LTISystem;
sys.Nx = 10; alpha = 0.8; rho = 1.5; actDens = 0.6; 
generate_dbl_stoch_chain(sys, rho, actDens, alpha);
sys.B1 = eye(sys.Nx);

tHorizon = 10;
x0       = rand(sys.Nx, 1);
w        = zeros(sys.Nx, tHorizon); % noiseless

% Params
params = MPCParams();
params.locality_ = 3;
params.tFIR_     = 3;
params.maxIters_ = 5000;
params.rho_      = 1; 
params.eps_p_    = 1e-3;
params.eps_d_    = 1e-3;
    
params.stateConsMtx_     = eye(sys.Nx);
params.stateUB_          = 1 * ones(sys.Nx, 1);
params.stateLB_          = -params.stateUB_;    
params.QSqrt_ = eye(sys.Nx);
params.RSqrt_ = eye(sys.Nu);

plotStates = [4 5];
plotInputs = [2 3];

%% Synthesize w/o terminal set
fprintf('Doing MPC with *no* terminal set...\n')
params.mode_        = MPCMode.Centralized;
[xCentA, uCentA, ~] = sls_mpc(sys, x0, w, params, tHorizon);

params.mode_        = MPCMode.Distributed;
[xA, uA, statsA]    = sls_mpc(sys, x0, w, params, tHorizon);

print_and_plot(params, xA, uA, xCentA, uCentA, 'Test A: No Terminal Set', plotStates, plotInputs);
fprintf('Test A: avgTime: %.4f, avgIters: %.4f\n\n', statsA.time_, statsA.iters_);

%% Add terminal set + synthesize w/ terminal set
fprintf('Synthesizing terminal set...\n')
[params, tStats] = terminal_set(sys, params);
params = remove_redundancy_terminal(sys, params);
fprintf('Terminal set: avgTime: %.4f, avgIters: %.4f\n\n', tStats.time_, tStats.iters_);

fprintf('Doing MPC *with* terminal set...\n')
params.mode_        = MPCMode.Centralized;
[xCentB, uCentB, ~] = sls_mpc(sys, x0, w, params, tHorizon);

params.mode_        = MPCMode.Distributed;
[xB, uB, statsB]    = sls_mpc(sys, x0, w, params, tHorizon);

print_and_plot(params, xB, uB, xCentB, uCentB, 'Test B: With Terminal Set', plotStates, plotInputs);
fprintf('Test B: avgTime: %.4f, avgIters: %.4f\n\n', statsB.time_, statsB.iters_);

