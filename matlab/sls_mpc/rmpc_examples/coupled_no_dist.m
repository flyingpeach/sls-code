clear all; close all; clc;

%% Setup plant + parameters
rng(420);
sys    = LTISystem;
sys.Nx = 4; Nx = sys.Nx; sys.B1 = eye(Nx);
alpha = 0.8; rho = 2; actDens = 0.7; 
generate_dbl_stoch_chain(sys, rho, actDens, alpha);

tHorizon = 10;
x0       = rand(sys.Nx, 1);
w        = zeros(sys.Nx, tHorizon); % noiseless

params = MPCParams();
params.locality_ = 3;
params.tFIR_     = 5;
params.maxIters_ = 5000;
params.rho_      = 1; 
params.eps_p_    = 5e-3;
params.eps_d_    = 5e-3;

% Coupling in both Q and R
params.QSqrt_ = diag(ones(Nx,1)) + diag(-1/3*ones(Nx-2,1),2) + diag(-1/2*ones(Nx-2,1),-2);
params.RSqrt_ = diag(randi([1, 3], sys.Nu, 1));
params.RSqrt_(1,2) = 1;

plotStates = 1:2;
plotInputs = 1:2;

params.distConsMtx_ = eye(sys.Nx);
params.distUB_      = zeros(sys.Nx, 1);
params.distLB_      = zeros(sys.Nx, 1);

%% Coupling in state and input constraints
% Adapted from tests_algorithm_2 in main toolbox
KState = zeros(sys.Nx);
KState(1,:) = [1 1 0 0];
KState(2,:) = [0 1 0 0];

params.stateConsMtx_ = KState;
params.stateUB_      = 2 * ones(sys.Nx, 1);
params.stateLB_      = inf(sys.Nx, 1);

KInput = zeros(sys.Nu);
KInput(1,:) = [1 1 0];

params.inputConsMtx_ = KInput;
params.inputUB_      = inf(sys.Nu, 1);
params.inputLB_      = -0.8* ones(sys.Nu, 1); % tight constraint

params.mode_        = MPCMode.Centralized;
[xCentA, uCentA, ~] = sls_mpc(sys, x0, w, params, tHorizon);

params.mode_        = MPCMode.Distributed;
[xA, uA, statsA]    = sls_mpc(sys, x0, w, params, tHorizon);

print_and_plot(params, xA, uA, xCentA, uCentA, 'Coupled rMPC', plotStates, plotInputs);
fprintf('avgTime: %.4f, avgIters: %.4f\n\n', statsA.time_, statsA.iters_);
