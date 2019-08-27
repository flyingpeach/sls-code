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
params       = SLSParams;
params.tFIR_ = 20;
params.obj_  = Objective.H2; % objective function

% simulation parameters
TMax                  = 25;                  % amount of time to simulate
w                     = zeros(sys.Nx, TMax); % disturbance
w(floor(sys.Nx/2), 1) = 10;

%% (1) basic sls with no constraints (centralized controller)
params.mode_      = SLSMode.Basic; % d-localized

[R1, M1, clnorm1] = state_fdbk_sls(sys, params);

[x1, u1] = simulate_system(sys, params, TMax, R1, M1, w);
plot_heat_map(x1, sys.B2*u1, 'Centralized');

%% (2) d-localized sls
params.mode_      = SLSMode.DLocalized;
params.actDelay_  = 1;
params.cSpeed_    = 2; % communication speed must be sufficiently large
params.d_         = 3;

[R2, M2, clnorm2] = state_fdbk_sls(sys, params);

[x2, u2] = simulate_system(sys, params, TMax, R2, M2, w);
plot_heat_map(x2, sys.B2*u2, 'Localized');

%% (3) approximate d-localized sls
params.mode_      = SLSMode.ApproxDLocalized;
params.cSpeed_    = 1;
params.lambda_    = 10^3;

[R3, M3, clnorm3, robust_stab] = state_fdbk_sls(sys, params);

[x3, u3] = simulate_system(sys, params, TMax, R3, M3, w);
plot_heat_map(x3, sys.B2*u3, 'Approximately Localized');