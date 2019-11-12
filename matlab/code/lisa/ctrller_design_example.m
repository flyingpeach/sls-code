%% solve for R, M first
acc_setup;

% for rank/zero conditions, try to match the precision of cvx_precision low
% http://cvxr.com/cvx/doc/solver.html#solver-precision
eps_rank    = 2.22e-16;
eps_nullsp_ = eps_rank.^(3/8);

%% sandbox
cParams = CtrllerParams();

% Basic, EncourageDelay, EncouargeLocal
% Delayed, Localized, DAndL
cParams.mode_ = SLSMode.Basic;

cParams.eps_nullsp_ = eps_nullsp_;
cParams.tFIR_       = slsParams.tFIR_;

% uncomment as needed
%cParams.actDelay_    = 0;
%cParams.cSpeed_      = 2;
%cParams.d_           = 4;
%cParams.CLDiffPen_   = 1e2;
%cParams.fastCommPen_ = 1e2;
%cParams.nonLocalPen_ = 1e2;

ctrller = find_ctrller(sys, slsParams, slsOuts, cParams);

%% visualize
visualize_RM_RMc(sys, slsOuts, ctrller, 'all');

%% calculate stats for new controller
cStats = CtrllerStats(eps_rank, cParams.tFIR_);
cs{1}  = ctrller;
cStats = get_ctrller_stats(cStats, sys, slsParams, slsOuts, cs);
cStats

%% parameter sweep (Tc)
cParams.mode_ = SLSMode.Basic;

Tcs      = 17:22;
numTcs   = length(Tcs);
csSweep = cell(numTcs, 1);
for idx=1:numTcs
    cParams.tFIR_ = Tcs(idx);
    csSweep{idx}  = find_ctrller(sys, slsParams, slsOuts, cParams);
end

%% parameter sweep (non-Tc)
cParams.tFIR_ = slsParams.tFIR_;
cParams.mode_ = SLSMode.Localized;

cParams.CLDiffPen_ = 1e2;

params    = 4:10;
numParams = length(params);
csSweep   = cell(numParams, 1);

for idx=1:numParams
    cParams.d_   = params(idx);
    csSweep{idx} = find_ctrller(sys, slsParams, slsOuts, cParams);
end

%% plotter for parameter sweep
% we might not want to plot all resultant ctrllers; can adjust below

% user input %%%%%%%%%%%%%%%%%%%%%%% 
sweepParamName = 'd-locality';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(sweepParamName, 'Tc')
    xSeries = Tcs; xSize = numTcs;
else
    xSeries = params; xSize = numParams;    
end

% user input %%%%%%%%%%%%%%%%%%%%%%% 
xWanted = xSeries;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

myIdx   = [];
for i=1:xSize
    x = xSeries(i);
    % ismember doesn't work well with floats    
    if find(abs(xWanted-x) < eps_rank) 
        myIdx = [myIdx, i];
    end
end

xPlot       = xSeries(myIdx);
csSweepPlot = csSweep(myIdx);

if strcmp(sweepParamName, 'Tc')
    cStats = CtrllerStats(eps_rank, xPlot);
else
    cStats = CtrllerStats(eps_rank, cParams.tFIR_, sweepParamName, xPlot);
end

cStats = get_ctrller_stats(cStats, sys, slsParams, slsOuts, csSweepPlot);

% plot metrics and save to file
savepath = 'C:\Users\flyin\Desktop\caltech\research\sls controller design';
plot_ctrller_stats(cStats, savepath);
