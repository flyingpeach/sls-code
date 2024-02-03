clear all; close all; clc;

%% Setup plant + parameters
rng(420);
Nx     = 8; alpha = 0.8; rho = 2; actDens = 0.6; 
sys    = generate_dbl_stoch_chain(Nx, rho, actDens, alpha);
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

params.QSqrt_ = diag(randi([1, 3], sys.Nx, 1));
params.RSqrt_ = diag(randi([1, 3], sys.Nu, 1));

plotStates = 1;
plotInputs = 2;

%% TEST A: Nominal, no coupling, no constraints
params.mode_        = MPCMode.Distributed;
[xA, uA, statsA]    = sls_mpc(sys, x0, w, params, tHorizon);

params.mode_        = MPCMode.Centralized;
[xCentA, uCentA, ~] = sls_mpc(sys, x0, w, params, tHorizon);

print_and_plot(params, xA, uA, xCentA, uCentA, 'Uncoupled Test A', plotStates, plotInputs);
fprintf('avgTime: %.4f, avgIters: %.4f\n\n', statsA.time_, statsA.iters_);

%% TEST B: Nominal, no coupling, state constraints
params.stateConsMtx_ = eye(sys.Nx);
params.stateUB_      =  1.7 * ones(sys.Nx, 1);  % not a tight constraint
params.stateLB_      = -0.5 * ones(sys.Nx, 1); % tight constraint

params.mode_        = MPCMode.Distributed;
[xB, uB, statsB]    = sls_mpc(sys, x0, w, params, tHorizon);

params.mode_        = MPCMode.Centralized;
[xCentB, uCentB, ~] = sls_mpc(sys, x0, w, params, tHorizon);

print_and_plot(params, xB, uB, xCentB, uCentB, 'Uncoupled Test B', plotStates, plotInputs);
fprintf('avgTime: %.4f, avgIters: %.4f\n\n', statsB.time_, statsB.iters_);

%% TEST C: Nominal, no coupling, state + input constraints
% state constraints still apply from TEST B if run sequentially
params.inputConsMtx_ = eye(sys.Nu);
params.inputUB_      = inf(sys.Nu, 1);
params.inputLB_      = -1.5 * ones(sys.Nu, 1); % tight constraint

params.mode_        = MPCMode.Distributed;
[xC, uC, statsC]    = sls_mpc(sys, x0, w, params, tHorizon);

params.mode_        = MPCMode.Centralized;
[xCentC, uCentC, ~] = sls_mpc(sys, x0, w, params, tHorizon);

print_and_plot(params, xC, uC, xCentC, uCentC, 'Uncoupled Test C', plotStates, plotInputs);
fprintf('avgTime: %.4f, avgIters: %.4f\n\n', statsC.time_, statsC.iters_);
