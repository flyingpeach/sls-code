%% solve for R, M
clear; close all; clc; 

% specify system matrices
sys    = LTISystem();
sys.Nx = 10;

alpha = 0.4; rho = 1; actDens = 0.3;
generate_dbl_stoch_chain(sys, rho, actDens, alpha); % generate sys.A, sys.B2

sys.B1  = eye(sys.Nx);
sys.C1  = [speye(sys.Nx); sparse(sys.Nu, sys.Nx)];
sys.D12 = [sparse(sys.Nx, sys.Nu); speye(sys.Nu)];

% sls parameters
slsParams           = SLSParams();
slsParams.T_        = 20;
slsParams.obj_      = Objective.H2; % objective function
slsParams.mode_     = SLSMode.Basic;
slsParams.robCoeff_ = 1e3;

slsOuts = state_fdbk_sls(sys, slsParams);

% for rank/zero conditions, try to match the precision of cvx_precision low
% http://cvxr.com/cvx/doc/solver.html#solver-precision
eps_base   = 2.22e-16;
eps_nullsp = eps_base.^(3/8);

%% sandbox
cParams = CtrllerParams();

% Basic, EncourageDelay, EncouargeLocal
% Delayed, Localized, DAndL
cParams.mode_ = SLSMode.Basic;

cParams.eps_nullsp_ = eps_nullsp;
cParams.T_          = slsParams.T_;

% uncomment as needed
%cParams.actDelay_    = 0;
%cParams.cSpeed_      = 2;
%cParams.d_           = 4;
%cParams.CLDiffPen_   = 1e2;
%cParams.fastCommPen_ = 1e2;
%cParams.nonLocalPen_ = 1e2;

ctrller = find_ctrller(sys, slsOuts, cParams);

%% visualize
plot_ctrller_matrices(sys, slsOuts, ctrller, 'all');

%% calculate stats for new controller
cStats = CtrllerStats(eps_nullsp);
cs{1}  = ctrller;
cStats = get_ctrller_stats(cStats, sys, slsOuts, cs);
cStats

%% parameter sweep
params    = 19:21;
numParams = length(params);
csSweep   = cell(numParams, 1);

for idx=1:numParams
    cParams.T_   = params(idx);
    csSweep{idx} = find_ctrller(sys, slsOuts, cParams);
end

%% plotter for parameter sweep
% we might not want to plot all resultant ctrllers; can adjust below

% user input %%%%%%%%%%%%%%%%%%%%%%% 
sweepParamName = 'Tc';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xSeries = params; 
xSize   = numParams;    

% user input %%%%%%%%%%%%%%%%%%%%%%% 
xWanted = xSeries;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

myIdx   = [];
for i=1:xSize
    x = xSeries(i);
    % ismember doesn't work well with floats    
    if find(abs(xWanted-x) < eps_base) 
        myIdx = [myIdx, i];
    end
end
xPlot       = xSeries(myIdx);
csSweepPlot = csSweep(myIdx);

cStats = CtrllerStats(eps_nullsp, sweepParamName, xPlot);
cStats = get_ctrller_stats(cStats, sys, slsOuts, csSweepPlot);

% plot metrics and save to file
savepath = 'C:\Users\flyin\Desktop\caltech\research\sls controller design';
plot_ctrller_stats(cStats, savepath);
