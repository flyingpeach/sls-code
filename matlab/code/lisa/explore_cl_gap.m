%% solve original SLS problem
clear; close all; clc; 

% set up system / parameters
sys    = LTISystem();
sys.Nx = 10;

alpha = 0.4; rho = 1; actDens = 0.3;
generate_dbl_stoch_chain(sys, rho, actDens, alpha); % generate sys.A, sys.B2

sys.B1  = eye(sys.Nx);
sys.C1  = [speye(sys.Nx); sparse(sys.Nu, sys.Nx)];
sys.D12 = [sparse(sys.Nx, sys.Nu); speye(sys.Nu)];

slsParams       = SLSParams();
slsParams.T_    = 20;
slsParams.obj_  = Objective.H2;

% centralized SLS
slsParams.mode_ = SLSMode.Basic;
slsOuts         = state_fdbk_sls(sys, slsParams);

%% parameter sweep
eps_base   = 2.22e-16;
eps_nullsp = eps_base.^(3/8);

cParams             = CtrllerParams();
cParams.mode_       = slsParams.mode_;
cParams.eps_nullsp_ = eps_nullsp;
cParams.T_          = slsParams.T_;

params    = linspace(eps_base.^(4/8), eps_base.^(0.1/8), 15);
numParams = length(params);
csSweep   = cell(numParams, 1);

gaps = zeros(numParams, 1);
for idx=1:numParams
    cParams.eps_nullsp_ = params(idx)/8;
    csSweep{idx} = find_ctrller(sys, slsOuts, cParams);

    gaps(idx)    = get_constr_gap(sys, cParams, slsOuts, csSweep{idx});
end

cStats = CtrllerStats(eps_nullsp, 'eps_null_relaxation', params);
cStats = get_ctrller_stats(cStats, sys, slsOuts, csSweep);

figure(1)
subplot(2,1,1)
hold on
plot(params, cStats.GcDiffs + cStats.HcDiffs);
plot(params, gaps);
legend('CLMapGap', 'ConstrGap');

subplot(2,1,2)
plot(gaps, cStats.GcDiffs + cStats.HcDiffs);
title('CLMapGap vs ConstrGap');