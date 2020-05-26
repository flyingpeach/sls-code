clear all; close all; clc;

%% Setup plant + parameters
rng(420);
sys    = LTISystem;
sys.Nx = 8; 
alpha = 0.8; rho = 2; actDens = 0.6; 
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

params.QSqrt_ = diag(randi([1, 3], sys.Nx, 1));
params.RSqrt_ = diag(randi([1, 3], sys.Nu, 1));

params.mode_  = MPCMode.Distributed;

plotStates = 1;
plotInputs = 2;

%% TEST A: Algorithm 1, no constraints
paramsCent       = copy(params);
paramsCent.mode_ = MPCMode.Centralized;

xA = zeros(sys.Nx, tHorizon); xCentA = zeros(sys.Nx, tHorizon);
uA = zeros(sys.Nu, tHorizon); uCentA = zeros(sys.Nu, tHorizon);
xA(:,1) = x0; xCentA(:,1) = x0;
for t=1:tHorizon-1
    fprintf('Calculating time %d of %d\n', t, tHorizon)
    [xA(:,t+1), uA(:,t), ~]           = sls_mpc(sys, xA(:,t), params);
    [xCentA(:,t+1), uCentA(:,t), ~] = sls_mpc(sys, xCentA(:,t), paramsCent);
end

print_and_plot(params, xA, uA, xCentA, uCentA, 'Alg1 Test A', plotStates, plotInputs);

%% TEST B: Algorithm 1, with state constraints
params.stateConsMtx_ = eye(sys.Nx);
params.stateUB_      = 1.7;  % not a tight constraint
params.stateLB_      = -0.5; % tight constraint

paramsCent       = copy(params);
paramsCent.mode_ = MPCMode.Centralized;

xB = zeros(sys.Nx, tHorizon); xCentB = zeros(sys.Nx, tHorizon);
uB = zeros(sys.Nu, tHorizon); uCentB = zeros(sys.Nu, tHorizon);
xB(:, 1) = x0; xCentB(:, 1) = x0;
for t=1:tHorizon-1
    fprintf('Calculating time %d of %d\n', t, tHorizon)
    [xB(:,t+1), uB(:,t), ~]           = sls_mpc(sys, xB(:,t), params);
    [xCentB(:,t+1), uCentB(:,t), ~] = sls_mpc(sys, xCentB(:,t), paramsCent);
end

print_and_plot(params, xB, uB, xCentB, uCentB, 'Alg1 Test B', plotStates, plotInputs);

%% TEST C: Algorithm 1, with state + input constraints
% state constraints still apply from TEST B if run sequentially
params.inputConsMtx_ = eye(sys.Nu);
params.inputLB_      = -1.5; % tight constraint

paramsCent       = copy(params);
paramsCent.mode_ = MPCMode.Centralized;

xC = zeros(sys.Nx, tHorizon); xCentC = zeros(sys.Nx, tHorizon);
uC = zeros(sys.Nu, tHorizon); uCentC = zeros(sys.Nu, tHorizon);
xC(:, 1) = x0; xCentC(:, 1) = x0;
for t=1:tHorizon-1
    fprintf('Calculating time %d of %d\n', t, tHorizon)
    [xC(:,t+1), uC(:,t), ~]           = sls_mpc(sys, xC(:,t), params);
    [xCentC(:,t+1), uCentC(:,t), ~] = sls_mpc(sys, xCentC(:,t), paramsCent);
end

print_and_plot(params, xC, uC, xCentC, uCentC, 'Alg1 Test C', plotStates, plotInputs);
