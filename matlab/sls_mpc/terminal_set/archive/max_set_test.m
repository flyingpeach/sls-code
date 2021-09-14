% Many called functions use the mpt3 toolbox! Please install it
% Copied from real simulations
clear; clc; close;

%% 
gridSize      = 3; % total size = gridSize^2
seed          = randi([0 100]);
connectThresh = 0.65; % lower = more connected grid

numNodes      = gridSize * gridSize;

% Actuation
actDens        = 0.8;
numActs        = round(actDens*numNodes);
actuatedNodes  = randsample(numNodes, numActs);

% Sampling time
Ts = 0.2; 

[adjMtx, nodeCoords, susceptMtx, inertiasInv, dampings] = setup_grid_topology(gridSize, connectThresh, seed);
sys = setup_grid_plant(actuatedNodes, adjMtx, susceptMtx, inertiasInv, dampings, Ts);

params = MPCParams();

% Input and state constraints
params.stateConsMtx_ = eye(sys.Nx);
params.stateUB_      = randi(8, sys.Nx, 1);
params.stateLB_      = -params.stateUB_;

params.inputConsMtx_ = eye(sys.Nu);
params.inputUB_      = randi(8, sys.Nu, 1);
params.inputLB_      = -params.inputUB_;

[Hmax_, hmax_, iters] = nominal_max_ci_set(sys, params);
iters