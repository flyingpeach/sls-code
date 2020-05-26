clear all; close all; clc;

%% Setup plant + parameters
rng(420);
sys2    = LTISystem;
sys2.Nx = 4; 
alpha = 0.8; rho = 2; actDens = 0.7; 
generate_dbl_stoch_chain(sys2, rho, actDens, alpha);

x0 = rand(sys2.Nx, 1);

params = MPCParams();
params.locality_ = 3;
params.tFIR_     = 5;
params.tHorizon_ = 10;
params.maxIters_ = 5000;
params.rho_      = 1; 
params.eps_p_    = 1e-3;
params.eps_d_    = 1e-3;

params.maxItersCons_ = 500;
params.mu_           = 1;
params.eps_x_        = 1e-3;
params.eps_z_        = 1e-3;

Nx = sys2.Nx;
params.QSqrt_ = diag(ones(Nx,1)) + diag(-1/2*ones(Nx-2,1),2) + diag(-1/2*ones(Nx-2,1),-2);
params.RSqrt_ = diag(randi([1, 3], sys2.Nu, 1));

% State constraints (different coupling from cost)
K = zeros(Nx);
K(1,:) = [1 1 0 0];
K(2,:) = [0 1 1 0];

plotStates = [2 3];
plotInputs = 1:sys2.Nu;

%% TEST A: Algorithm 2, no constraints
params.mode_    = MPCMode.Distributed;
[x, u, ~]       = sls_mpc(sys2, x0, params);

params.mode_    = MPCMode.Centralized;
[xVal, uVal, ~] = sls_mpc(sys2, x0, params);

print_and_plot(params, x, u, xVal, uVal, 'Alg2 Test A', plotStates, plotInputs);

%% TEST B: Algorithm 2, with state constraints
params.stateConsMtx_ = K;
params.stateUB_      = 0.8; % tight constraint

params.mode_    = MPCMode.Distributed;
[x, u, ~]       = sls_mpc(sys2, x0, params);

params.mode_    = MPCMode.Centralized;
[xVal, uVal, ~] = sls_mpc(sys2, x0, params);

print_and_plot(params, x, u, xVal, uVal, 'Alg2 Test B', plotStates, plotInputs);

%% TEST C: Algorithm 2, with state + input constraints
% state constraints still apply from TEST B if run sequentially
params.inputConsMtx_ = eye(sys2.Nu);
params.inputLB_      = -0.8; % tight constraint

params.mode_   = MPCMode.Distributed;
[x, u, ~]       = sls_mpc(sys2, x0, params);

params.mode_   = MPCMode.Centralized;
[xVal, uVal, ~] = sls_mpc(sys2, x0, params);

print_and_plot(params, x, u, xVal, uVal, 'Alg2 Test C', plotStates, plotInputs);
