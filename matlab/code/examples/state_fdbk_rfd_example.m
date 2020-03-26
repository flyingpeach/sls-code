% note: run sequentially; some params from later sections depend on
% params specified in previous sections
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
slsParams      = SLSParams();
slsParams.T_   = 15;
slsParams.obj_ = Objective.H2; % objective function

% delayed and localized sls with rfd
num_acts = []; clnorms = [];

slsParams.mode_      = SLSMode.DAndL;
slsParams.actDelay_  = 1;
slsParams.cSpeed_    = 2;
slsParams.d_         = 3;

for power = -2:1:3
    slsParams.rfdCoeff_ = 10^power;
    slsParams.rfd_      = true;
    slsOutsRFD2         = state_fdbk_sls(sys, slsParams);

     % check performance with rfd-designed system
    sysAfterRFD2     = updateActuation(sys, slsOutsRFD2);
    slsParams.rfd_   = false;
    slsOutsAfterRFD2 = state_fdbk_sls(sysAfterRFD2, slsParams);

    num_acts         = [num_acts; length(slsOutsRFD2.acts_)];
    clnorms          = [clnorms; slsOutsAfterRFD2.clnorm_];
end

figure;
plot(num_acts,clnorms,'*-');
xlabel('Number of actuators'); ylabel('Close loop norm');
title('d-localized RFD tradeoff curve');
