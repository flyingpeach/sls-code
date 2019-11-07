%% choose the system you want to work with
setup1;

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
cParams.tc_         = slsParams.tFIR_;

% uncomment as needed
% cParams.actDelay_    = 0;
% cParams.cSpeed_      = 0;
% cParams.d_           = 0;
% cParams.CLDiffPen_   = 0;
% cParams.fastCommPen_ = 0;
% cParams.nonLocalPen_ = 0;

ctrller = find_ctrller(sys, slsParams, slsOuts, cParams);

%% visualize
visualize_RM_RMc(sys, slsOuts, ctrller, 'all');




%% calculate metrics
zThresh = tol;
met     = AltImplMetrics(zThresh, Tc);
met     = calc_mtx_metrics(met, sys, slsParams, slsOuts, s_a);
met     = calc_cl_metrics(met, sys, simParams, slsParams, slsOuts, s_a);
met

%% parameter sweep (Tc)
Tcs    = 2:18;
numTcs = length(Tcs);
slsOuts_alts = cell(numTcs, 1);
for idx=1:numTcs
    Tc = Tcs(idx);
    slsOuts_alts{idx} = find_alt_impl(sys, slsParams, slsOuts, Tc, settings);
end

%% parameter sweep (others)
Tc = round(slsParams.tFIR_/2);

settings.mode_ = AltImplMode.ApproxLS;
params       = [1/12, 1/8, 1/6, 1/4];
numParams    = length(params);
slsOuts_alts = cell(numParams, 1);

for i=1:numParams
    settings.svThresh_ = eps.^(params(i));
    slsOuts_alts{i}  = find_alt_impl(sys, slsParams, slsOuts, Tc, settings);
end

%% plot stuff
% we might not want to plot all of slsOuts_alts, so this is the sandbox to
% adjust which slsOuts_alts to plot

zThresh = tol;
sweepParamName = 'Tc';

if strcmp(sweepParamName, 'Tc')
    xSeries = Tcs;
    xSize   = numTcs;
else
    xSeries = params;
    xSize   = numParams;    
end

% can specify which x to plot
xWanted = xSeries;
myIdx   = [];

for i=1:xSize
    x = xSeries(i);
    if find(abs(xWanted-x)<eps) % ismember doesn't work well with floats
        myIdx = [myIdx, i];
    end
end

xPlot            = xSeries(myIdx);
slsOuts_altsPlot = slsOuts_alts(myIdx);

if strcmp(sweepParamName, 'Tc')
    met = AltImplMetrics(zThresh, xPlot);
else
    met = AltImplMetrics(zThresh, Tc, sweepParamName, xPlot);
end

% calculate matrix-specific metrics
met = calc_mtx_metrics(met, sys, slsParams, slsOuts, slsOuts_altsPlot);

% calculate closed-loop metrics and plot select heat maps
met = calc_cl_metrics(met, sys, simParams, slsParams, slsOuts, slsOuts_altsPlot);

% plot metrics and save to file
savepath = 'C:\Users\Lisa\Desktop\caltech\research\implspace\tmp\';
plot_metrics(met, savepath);

disp(['Statuses:', print_statuses(xSeries, slsOuts_alts)]); 

