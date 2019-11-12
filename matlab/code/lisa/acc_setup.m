% setup
clear; close all; clc; 

% specify system matrices
sys    = LTISystem;
sys.Nx = 10;

alpha = 0.4; rho = 1; actDens = 0.3;
% generate sys.A, sys.B2
generate_dbl_stoch_chain(sys, rho, actDens, alpha);

sys.B1  = eye(sys.Nx);
sys.C1  = [speye(sys.Nx); sparse(sys.Nu, sys.Nx)];
sys.D12 = [sparse(sys.Nx, sys.Nu); speye(sys.Nu)];

% sls parameters
slsParams           = SLSParams();
slsParams.tFIR_     = 20;
slsParams.obj_      = Objective.H2; % objective function
slsParams.mode_     = SLSMode.Basic;
slsParams.robCoeff_ = 1e3;

slsOuts = state_fdbk_sls(sys, slsParams);
