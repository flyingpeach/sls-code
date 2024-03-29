clear all; close all; clc;

%% Setup plant + parameters
rng(420);

Nx  = 4; alpha = 0.8; rho = 2; actDens = 0.7; 
sys = generate_dbl_stoch_chain(Nx, rho, actDens, alpha);
sys.B1 = eye(sys.Nx);

tHorizon = 10;
x0       = rand(sys.Nx, 1);
w        = zeros(sys.Nx, tHorizon); % noiseless

params = MPCParams();
params.locality_ = 3;
params.tFIR_     = 5;
params.maxIters_ = 5000;
params.rho_      = 1; 
params.eps_p_    = 1e-3;
params.eps_d_    = 1e-3;

params.maxItersCons_ = 500;
params.mu_           = 1;
params.eps_x_        = 1e-3;
params.eps_z_        = 1e-3;

Nx = sys.Nx;
params.QSqrt_ = diag(ones(Nx,1)) + diag(-1/2*ones(Nx-2,1),2) + diag(-1/2*ones(Nx-2,1),-2);
params.RSqrt_ = diag(randi([1, 3], sys.Nu, 1));

params.mode_  = MPCMode.Distributed;

% State constraints (different coupling from cost)
K = zeros(Nx);
K(1,:) = [1 1 0 0];
K(2,:) = [0 1 1 0];

plotStates = [2 3];
plotInputs = 1:sys.Nu;

%% TEST A: Nominal, coupled objective, no constraints
params.mode_        = MPCMode.Distributed;
[xA, uA, statsA]    = sls_mpc(sys, x0, w, params, tHorizon);

params.mode_        = MPCMode.Centralized;
[xCentA, uCentA, ~] = sls_mpc(sys, x0, w, params, tHorizon);

print_and_plot(params, xA, uA, xCentA, uCentA, 'Coupled Test A', plotStates, plotInputs);
fprintf('avgTime: %.4f, avgIters: %.4f\n\n', statsA.time_, statsA.iters_);

%% TEST B: Nominal, coupled objective + state constraints
params.stateConsMtx_ = K;
params.stateUB_      = 0.8 * ones(sys.Nx, 1); % tight constraint
params.stateLB_      = inf(sys.Nx, 1);

params.mode_        = MPCMode.Distributed;
[xB, uB, statsB]    = sls_mpc(sys, x0, w, params, tHorizon);

params.mode_        = MPCMode.Centralized;
[xCentB, uCentB, ~] = sls_mpc(sys, x0, w, params, tHorizon);

print_and_plot(params, xB, uB, xCentB, uCentB, 'Coupled Test B', plotStates, plotInputs);
fprintf('avgTime: %.4f, avgIters: %.4f\n\n', statsB.time_, statsB.iters_);

%% TEST C: Nominal, coupled objective + state constraints, input constraints
% state constraints still apply from TEST B if run sequentially
params.inputConsMtx_ = eye(sys.Nu);
params.inputUB_      = inf(sys.Nu, 1);
params.inputLB_      = -0.8 * ones(sys.Nu, 1); % tight constraint

params.mode_        = MPCMode.Distributed;
[xC, uC, statsC]    = sls_mpc(sys, x0, w, params, tHorizon);

params.mode_        = MPCMode.Centralized;
[xCentC, uCentC, ~] = sls_mpc(sys, x0, w, params, tHorizon);

print_and_plot(params, xC, uC, xCentC, uCentC, 'Coupled Test C', plotStates, plotInputs);
fprintf('avgTime: %.4f, avgIters: %.4f\n\n', statsC.time_, statsC.iters_);
