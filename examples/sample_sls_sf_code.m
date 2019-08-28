clear; close all; clc; 

% specify system matrices
sys    = LTISystem;
sys.Nx = 10;

alpha = 0.2; rho = 1; actDens = 1;
generate_dbl_stoch_chain(sys, rho, actDens, alpha); % generate sys.A, sys.B2

sys.B1  = eye(sys.Nx); % used in simulation
sys.C1  = [speye(sys.Nx); sparse(sys.Nu, sys.Nx)]; % used in H2/HInf ctrl
sys.D12 = [sparse(sys.Nx, sys.Nu); speye(sys.Nu)];

% sls parameters
slsParams       = SLSParams;
slsParams.tFIR_ = 20;
slsParams.obj_  = Objective.H2; % objective function

% simulation parameters
TMax                  = 25;                  % amount of time to simulate
w                     = zeros(sys.Nx, TMax); % disturbance
w(floor(sys.Nx/2), 1) = 10;

%% (1) basic sls (centralized controller)
slsParams.mode_ = SLSMode.Basic;

slsOuts1 = state_fdbk_sls(sys, slsParams);
[x1, u1] = simulate_system(sys, slsParams, slsOuts1, TMax, w);
plot_heat_map(x1, sys.B2*u1, 'Centralized');

%% (2) d-localized sls
slsParams.mode_     = SLSMode.DLocalized;
slsParams.actDelay_ = 1;
slsParams.cSpeed_   = 2; % communication speed must be sufficiently large
slsParams.d_        = 3;

slsOuts2 = state_fdbk_sls(sys, slsParams);
[x2, u2] = simulate_system(sys, slsParams, slsOuts2, TMax, w);
plot_heat_map(x2, sys.B2*u2, 'Localized');

%% (3) approximate d-localized sls
slsParams.mode_     = SLSMode.ApproxDLocalized;
slsParams.cSpeed_   = 1;
slsParams.robCoeff_ = 10^3;

slsOuts3 = state_fdbk_sls(sys, slsParams);
[x3, u3] = simulate_system(sys, slsParams, slsOuts3, TMax, w);
plot_heat_map(x3, sys.B2*u3, 'Approximately Localized');