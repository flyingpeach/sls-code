clear; close all; clc;

% specify system matrices
sys    = LTISystem;
sys.Nx = 10;

alpha = 0.2; rho = 0.8; actDens = 1;
randn('seed', 0);
generate_rand_chain(sys, rho, actDens); % generate sys.A, sys.B2

sys.B1  = eye(sys.Nx); % used in simulation
sys.C1  = [speye(sys.Nx); sparse(sys.Nu, sys.Nx)]; % used in H2/HInf ctrl
sys.D12 = [sparse(sys.Nx, sys.Nu); speye(sys.Nu)];

% sls parameters
params       = SLSParams;
params.tFIR_ = 15;
params.obj_  = Objective.H2; % objective function

% simulation parameters
TMax                  = 25;                  % amount of time to simulate
w                     = zeros(sys.Nx, TMax); % disturbance
w(floor(sys.Nx/2), 1) = 10;

%% (1) basic sls (centralized controller) with rfd
num_acts = []; clnorms = [];

params.mode_ = SLSMode.Basic;
sysAfterRFD  = copy(sys); % contains actuation matrices designed by rfd
   
for power = -2:1:3
    params.rfdCoeff_ = 10^power;
    [R1, M1, acts]   = state_fdbk_sls_rfd(sys, params);   
    sysAfterRFD.B2   = sys.B2(:, acts);
    sysAfterRFD.D12  = sys.D12(:, acts);
    sysAfterRFD.Nu   = size(acts, 1);
    [R1, M1, clnorm] = state_fdbk_sls(sysAfterRFD, params); % find clnorm 
    num_acts         = [num_acts; length(acts)];
    clnorms          = [clnorms; clnorm];
end

figure
plot(num_acts,clnorms,'*-');
xlabel('Number of actuators')
ylabel('Close loop norm')
title('Centralized RFD tradeoff curve')

%% (2) d-localized sls with rfd
num_acts = []; clnorms = [];

params.mode_      = SLSMode.DLocalized;
params.actDelay_  = 1;
params.cSpeed_    = 2;
params.d_         = 3;

for power = -2:1:3
    params.rfdCoeff_ = 10^power;
    [R2, M2, acts]   = state_fdbk_sls_rfd(sys, params);   
    sysAfterRFD.B2   = sys.B2(:, acts);
    sysAfterRFD.D12  = sys.D12(:, acts);
    sysAfterRFD.Nu   = size(acts, 1);
    [R2, M2, clnorm] = state_fdbk_sls(sysAfterRFD, params); % find clnorm 
    num_acts         = [num_acts; length(acts)];
    clnorms          = [clnorms; clnorm];
end

figure
plot(num_acts,clnorms,'*-');
xlabel('Number of actuators')
ylabel('Close loop norm')
title('d-localized RFD tradeoff curve')

%% (3) approximate d-localized sls with rfd
num_acts = []; clnorms = [];

params.mode_      = SLSMode.ApproxDLocalized;
params.robCoeff_  = 10^4;

for power = -2:1:3
    params.rfdCoeff_ = 10^power;
    [R3, M3, acts]   = state_fdbk_sls_rfd(sys, params);
    sysAfterRFD.B2   = sys.B2(:, acts);
    sysAfterRFD.D12  = sys.D12(:, acts);
    sysAfterRFD.Nu   = size(acts, 1);
    [R3, M3, clnorm] = state_fdbk_sls(sysAfterRFD, params); % find clnorm 
    num_acts         = [num_acts; length(acts)];
    clnorms          = [clnorms; clnorm];
end

figure
plot(num_acts,clnorms,'*-');
xlabel('Number of actuators')
ylabel('Close loop norm')
title('Approx d-localized RFD tradeoff curve')