clear; close all; clc; 

% which things we want to plot
plotAnimation = false;
plotHeatMap   = true;

% specify system matrices
sys    = LTISystem;
sys.Nx = 50;

alpha = 0.8; rho = 1; actDens = 1; 
generate_dbl_stoch_chain(sys, rho, actDens, alpha); % update A, B2
sys.B1  = eye(sys.Nx); % used in simulation
sys.C1  = [speye(sys.Nx); sparse(sys.Nu, sys.Nx)];
sys.D12 = [sparse(sys.Nx, sys.Nu); speye(sys.Nu)];

% sls parameters
slsParams           = SLSParams();
slsParams.obj_      = Objective.H2;
slsParams.T_        = 10;
slsParams.actDelay_ = 1;
slsParams.cSpeed_   = 1;
slsParams.d_        = 8;
slsParams.robCoeff_ = 1000;
slsParams.mode_     = SLSMode.DAndL;
slsParams.approx_   = true;

% simulation parameters
simParams           = SimParams;
simParams.tSim_     = 80;
simParams.openLoop_ = false;
simParams.w_        = zeros(sys.Nx, 100);
simParams.w_(floor(sys.Nx/2), 1) = 10;

slsOuts = state_fdbk_sls(sys, slsParams);
[x, u]  = simulate_system(sys, simParams, slsOuts.R_, slsOuts.M_);

nodeCoords = [1:1:sys.Nx;
              zeros(1, sys.Nx)]';

onesVec = ones(sys.Nx - 1, 1);
adjMtx  = eye(sys.Nx) + diag(onesVec, 1) + diag(onesVec, -1);

if plotAnimation
    waitTime = 0.01;
    logScale = true; % want colours to match the heat map colours
    plot_graph_animation(adjMtx, nodeCoords, slsParams, x, sys.B2*u, waitTime, logScale);
end

if plotHeatMap
    plot_heat_map(x, sys.B2*u, '');
end