% setup1 (basic / centralized)
%clear; close all; clc; 

% specify system matrices
sys    = LTISystem;
sys.Nx = 10;

alpha = 0.4; rho = 1; actDens = 0.3;
generate_dbl_stoch_chain(sys, rho, actDens, alpha); % generate sys.A, sys.B2

sys.B1  = eye(sys.Nx); % used in simulation
sys.C1  = [speye(sys.Nx); sparse(sys.Nu, sys.Nx)]; % used in H2/HInf ctrl
sys.D12 = [sparse(sys.Nx, sys.Nu); speye(sys.Nu)];

% sls parameters
slsParams       = SLSParams;
slsParams.tFIR_ = 20;
slsParams.obj_  = Objective.H2; % objective function
slsParams.mode_ = SLSMode.ApproxDLocalized;
slsParams.robCoeff_ = 1e3;
%%
% simulation parameters
simParams           = SimParams;
simParams.tSim_     = 40;
simParams.w_        = zeros(sys.Nx, simParams.tSim_); % disturbance
simParams.w_(5, 1) = 1;
simParams.openLoop_ = false;

slsOuts = state_fdbk_sls(sys, slsParams);
[xOld, uOld] = simulate_system(sys, slsParams, slsOuts, simParams);
plot_heat_map(xOld, sys.B2*uOld, 'Original');
