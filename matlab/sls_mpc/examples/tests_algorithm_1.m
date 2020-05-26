clear all; close all; clc;

%% Setup plant + parameters
rng(420);
sys    = LTISystem;
sys.Nx = 8; 
alpha = 0.8; rho = 2; actDens = 0.6; 
generate_dbl_stoch_chain(sys, rho, actDens, alpha);

x0 = rand(sys.Nx, 1);

params = MPCParams();
params.locality_ = 3;
params.tFIR_     = 5;
params.tHorizon_ = 10;
params.maxIters_ = 5000;
params.rho_      = 1; 
params.eps_p_    = 1e-3;
params.eps_d_    = 1e-3;

params.QSqrt_ = diag(randi([1, 3], sys.Nx, 1));
params.RSqrt_ = diag(randi([1, 3], sys.Nu, 1));

plotStates = 1;
plotInputs = 2;

%% TEST A: Algorithm 1, no constraints
params.mode_    = MPCMode.Distributed;
[x, u, ~]       = sls_mpc(sys, x0, params);

params.mode_    = MPCMode.Centralized;
[xVal, uVal, ~] = sls_mpc(sys, x0, params);

print_and_plot(params, x, u, xVal, uVal, 'Alg1 Test A', plotStates, plotInputs);

%% TEST B: Algorithm 1, with state constraints
params.stateConsMtx_ = eye(sys.Nx);
params.stateUB_      = 1.7;  % not a tight constraint
params.stateLB_      = -0.5; % tight constraint

params.mode_    = MPCMode.Distributed;
[x, u, ~]       = sls_mpc(sys, x0, params);

params.mode_    = MPCMode.Centralized;
[xVal, uVal, ~] = sls_mpc(sys, x0, params);

print_and_plot(params, x, u, xVal, uVal, 'Alg1 Test B', plotStates, plotInputs);

%% TEST C: Algorithm 1, with state + input constraints
% state constraints still apply from TEST B if run sequentially
params.inputConsMtx_ = eye(sys.Nu);
params.inputLB_      = -1.5; % tight constraint

params.mode_    = MPCMode.Distributed;
[x, u, ~]       = sls_mpc(sys, x0, params);

params.mode_    = MPCMode.Centralized;
[xVal, uVal, ~] = sls_mpc(sys, x0, params);

print_and_plot(params, x, u, xVal, uVal, 'Alg1 Test C', plotStates, plotInputs);
