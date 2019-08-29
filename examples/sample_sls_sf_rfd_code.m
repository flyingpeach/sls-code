clear; close all; clc;

% specify system matrices
sys    = LTISystem;
sys.Nx = 10;

rho = 0.8; actDens = 1;
randn('seed', 0);
generate_rand_chain(sys, rho, actDens); % generate sys.A, sys.B2

sys.B1  = eye(sys.Nx); % used in simulation
sys.C1  = [speye(sys.Nx); sparse(sys.Nu, sys.Nx)]; % used in H2/HInf ctrl
sys.D12 = [sparse(sys.Nx, sys.Nu); speye(sys.Nu)];

% sls parameters
slsParams       = SLSParams;
slsParams.tFIR_ = 15;
slsParams.obj_  = Objective.H2; % objective function

%% (1) basic sls (centralized controller) with rfd
num_acts = []; clnorms = [];

slsParams.mode_ = SLSMode.Basic;
sysAfterRFD     = copy(sys); % contains actuation matrices designed by rfd
   
for power = -2:1:3
    slsParams.rfdCoeff_ = 10^power;
    slsParams.rfd_      = true;
    slsOutsRFD1         = state_fdbk_sls(sys, slsParams);
    
    % check performance with rfd-designed system
    sysAfterRFD.B2   = sys.B2(:, slsOutsRFD1.acts_);
    sysAfterRFD.D12  = sys.D12(:, slsOutsRFD1.acts_);
    sysAfterRFD.Nu   = size(slsOutsRFD1.acts_, 1);
    slsParams.rfd_   = false;
    slsOutsAfterRFD1 = state_fdbk_sls(sysAfterRFD, slsParams);

    num_acts         = [num_acts; length(slsOutsRFD1.acts_)];
    clnorms          = [clnorms; slsOutsAfterRFD1.clnorm_];
end

figure
plot(num_acts,clnorms,'*-');
xlabel('Number of actuators'); ylabel('Close loop norm');
title('Centralized RFD tradeoff curve');

%% (2) d-localized sls with rfd
num_acts = []; clnorms = [];

slsParams.mode_      = SLSMode.DLocalized;
slsParams.actDelay_  = 1;
slsParams.cSpeed_    = 2;
slsParams.d_         = 3;

for power = -2:1:3
    slsParams.rfdCoeff_ = 10^power;
    slsParams.rfd_      = true;
    slsOutsRFD2         = state_fdbk_sls(sys, slsParams);

     % check performance with rfd-designed system
    sysAfterRFD.B2   = sys.B2(:, slsOutsRFD2.acts_);
    sysAfterRFD.D12  = sys.D12(:, slsOutsRFD2.acts_);
    sysAfterRFD.Nu   = size(slsOutsRFD2.acts_, 1);
    slsParams.rfd_   = false;
    slsOutsAfterRFD2 = state_fdbk_sls(sysAfterRFD, slsParams);
   
    num_acts         = [num_acts; length(slsOutsRFD2.acts_)];
    clnorms          = [clnorms; slsOutsAfterRFD2.clnorm_];
end

figure;
plot(num_acts,clnorms,'*-');
xlabel('Number of actuators'); ylabel('Close loop norm');
title('d-localized RFD tradeoff curve');

%% (3) approximate d-localized sls with rfd
num_acts = []; clnorms = [];

slsParams.mode_      = SLSMode.ApproxDLocalized;
slsParams.robCoeff_  = 10^4;

for power = -2:1:3
    slsParams.rfdCoeff_ = 10^power;
    slsParams.rfd_      = true;
    slsOutsRFD3         = state_fdbk_sls(sys, slsParams);
    
    % check performance with rfd-designed system
    sysAfterRFD.B2   = sys.B2(:, slsOutsRFD3.acts_);
    sysAfterRFD.D12  = sys.D12(:, slsOutsRFD3.acts_);
    sysAfterRFD.Nu   = size(slsOutsRFD3.acts_, 1);
    slsParams.rfd_   = false;
    slsOutsAfterRFD3 = state_fdbk_sls(sysAfterRFD, slsParams);

    num_acts         = [num_acts; length(slsOutsRFD3.acts_)];
    clnorms          = [clnorms; slsOutsAfterRFD3.clnorm_];
end

figure;
plot(num_acts,clnorms,'*-');
xlabel('Number of actuators'); ylabel('Close loop norm');
title('Approx d-localized RFD tradeoff curve');