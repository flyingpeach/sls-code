clear all; close all; clc;

%% Setup plant + parameters
rng(420);
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
params.eps_p_    = 1e-3;
params.eps_d_    = 1e-3;

params.QSqrt_ = diag(randi([1, 3], sys.Nx, 1));
params.RSqrt_ = diag(randi([1, 3], sys.Nu, 1));

% Bounded disturbance
params.distConsMtx_ = eye(sys.Nx);
params.distUB_      = zeros(sys.Nx, 1);
params.distLB_      = zeros(sys.Nx, 1);

plotStates = 1;
plotInputs = 2;

%% State and input constraints
params.stateConsMtx_ = eye(sys.Nx);
params.stateUB_      =  1.7 * ones(sys.Nx, 1); % not a tight constraint
params.stateLB_      = -0.5 * ones(sys.Nx, 1); % tight constraint

params.inputConsMtx_ = eye(sys.Nu);
params.inputUB_      = inf(sys.Nu, 1);
params.inputLB_      = -1.5 * ones(sys.Nu, 1); % tight constraint

params.mode_        = MPCMode.Centralized;
[xCentA, uCentA, ~] = sls_mpc(sys, x0, w, params, tHorizon);

params.mode_     = MPCMode.Distributed;
[xA, uA, statsA] = sls_mpc(sys, x0, w, params, tHorizon);

print_and_plot(params, xA, uA, xCentA, uCentA, 'rMPC PolyNoise Test', plotStates, plotInputs);
fprintf('avgTime: %.4f, avgIters: %.4f\n\n', statsA.time_, statsA.iters_);