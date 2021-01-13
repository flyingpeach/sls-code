clear all; close all; clc;

%% Setup plant + parameters
rng(420);
sys    = LTISystem;
sys.Nx = 4; 
alpha = 0.8; rho = 2; actDens = 0.7; 
generate_dbl_stoch_chain(sys, rho, actDens, alpha);

tHorizon = 10;
x0       = rand(sys.Nx, 1);

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

%% TEST A: Algorithm 2, no constraints
params.mode_                  = MPCMode.Distributed;
[xA, uA, avgTimeA, avgItersA] = sls_mpc(sys, x0, params, tHorizon);

params.mode_        = MPCMode.Centralized;
[xCentA, uCentA, ~] = sls_mpc(sys, x0, params, tHorizon);

print_and_plot(params, xA, uA, xCentA, uCentA, 'Alg1 Test A', plotStates, plotInputs);
fprintf('\navgTime: %.4f, avgIters: %.4f\n', avgTimeA, avgItersA);

%% TEST B: Algorithm 2, with state constraints
params.stateConsMtx_ = K;
params.stateUB_      = 0.8 * ones(sys.Nx, 1); % tight constraint
params.stateLB_      = inf(sys.Nx, 1);

params.mode_                  = MPCMode.Distributed;
[xB, uB, avgTimeB, avgItersB] = sls_mpc(sys, x0, params, tHorizon);

params.mode_        = MPCMode.Centralized;
[xCentB, uCentB, ~] = sls_mpc(sys, x0, params, tHorizon);

print_and_plot(params, xB, uB, xCentB, uCentB, 'Alg1 Test B', plotStates, plotInputs);
fprintf('\navgTime: %.4f, avgIters: %.4f\n', avgTimeB, avgItersB);

%% TEST C: Algorithm 2, with state + input constraints
% state constraints still apply from TEST B if run sequentially
params.inputConsMtx_ = eye(sys.Nu);
params.inputUB_      = inf(sys.Nu, 1);
params.inputLB_      = -0.8 * ones(sys.Nu, 1); % tight constraint

params.mode_                  = MPCMode.Distributed;
[xC, uC, avgTimeC, avgItersC] = sls_mpc(sys, x0, params, tHorizon);

params.mode_        = MPCMode.Centralized;
[xCentC, uCentC, ~] = sls_mpc(sys, x0, params, tHorizon);

print_and_plot(params, xC, uC, xCentC, uCentC, 'Alg1 Test C', plotStates, plotInputs);
fprintf('\navgTime: %.4f, avgIters: %.4f\n', avgTimeC, avgItersC);
