clear all; close all; clc;

%% Setup plant + parameters
rng(425);
sys    = LTISystem;
sys.Nx = 8; sys.B1 = eye(sys.Nx);
alpha = 0.8; rho = 2; actDens = 0.6; 
generate_dbl_stoch_chain(sys, rho, actDens, alpha);

tHorizon = 10;
x0       = rand(sys.Nx, 1);
w        = zeros(sys.Nx, tHorizon);

params = MPCParams();
params.locality_ = 3;
params.tFIR_     = 5;
params.maxIters_ = 5000;
params.rho_      = 1; 
params.eps_p_    = 5e-3;
params.eps_d_    = 5e-3;

params.QSqrt_ = diag(randi([1, 3], sys.Nx, 1));
params.RSqrt_ = diag(randi([1, 3], sys.Nu, 1));

% Locally bounded disturbance
params.locNoiseBound_ = 1;

% Adaptive ADMM
params.tau_i_   = 2;
params.tau_d_   = 2;
params.muAdapt_ = 10;
params.rhoMax_  = 3;

plotStates = 3;
plotInputs = 2;

%% State constraints
params.stateConsMtx_ = eye(sys.Nx);
params.stateUB_      = inf(sys.Nx, 1);
params.stateLB_      = -inf(sys.Nx, 1);
params.stateLB_(3)   = -0.4; % without this, state < -0.4

params.mode_        = MPCMode.Centralized;
[xCentA, uCentA, ~] = sls_mpc(sys, x0, w, params, tHorizon);

params.mode_     = MPCMode.Distributed;
[xA, uA, statsA] = sls_mpc(sys, x0, w, params, tHorizon);

print_and_plot(params, xA, uA, xCentA, uCentA, 'Locally Bnd Noise, Test A', plotStates, plotInputs);
fprintf('avgTime: %.4f, avgIters: %.4f\n\n', statsA.time_, statsA.iters_);

%% Input constraints
% state constraints still apply from TEST B if run sequentially
params.inputConsMtx_ = eye(sys.Nu);
params.inputUB_      = 1.5*ones(sys.Nu, 1); % tight constraint
params.inputLB_      = -inf(sys.Nu, 1);

params.mode_        = MPCMode.Centralized;
[xCentB, uCentB, ~] = sls_mpc(sys, x0, w, params, tHorizon);

params.mode_     = MPCMode.Distributed;
[xB, uB, statsB] = sls_mpc(sys, x0, w, params, tHorizon);

print_and_plot(params, xB, uB, xCentB, uCentB, 'Locally Bnd Noise, Test B', plotStates, plotInputs);
fprintf('avgTime: %.4f, avgIters: %.4f\n\n', statsB.stime_, statsB.iters_);
