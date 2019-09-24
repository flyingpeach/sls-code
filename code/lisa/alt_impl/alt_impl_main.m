%% choose the system you want to work with
setup1;

% for rank/zero conditions, try to match the precision of cvx_precision low
% http://cvxr.com/cvx/doc/solver.html#solver-precision
eps = 2.22e-16;
tol = eps.^(3/8);
%close all;

settings = AltImplSettings();

%% sandbox
% ExactOpt, Analytic, ApproxLS, ApproxLeaky, StrictDelay, StrictLocal, EncourageDelay
settings.mode_ = AltImplMode.StrictLocal;

% uncomment as needed
settings.tol_          = eps.^(3/8);
settings.locality_     = 5;
settings.nonzeroPen_   = 1e0;
% settings.svThresh_    = eps.^(1/6);
% settings.delay_       = 1;
% settings.clDiffPen_   = 1e3;
% settings.fastCommPen_ = 1e0;

Tc = slsParams.tFIR_;

slsOuts_alt = find_alt_impl(sys, slsParams, slsOuts, Tc, settings);
s_a{1}      = slsOuts_alt;

%% get baseline comparison

[RZeros, MZeros] = get_locality_constraints(sys, settings.locality_);
slsOutsBaseline = copy(slsOuts);

for k=1:slsParams.tFIR_
    slsOutsBaseline.R_{k}(RZeros) = 0;
    slsOutsBaseline.M_{k}(MZeros) = 0;
end
s_a{1} = slsOutsBaseline;

%% calculate metrics
zThresh = tol;
met     = AltImplMetrics(zThresh, Tc);
met     = calc_mtx_metrics(met, sys, slsParams, slsOuts, s_a);
met     = calc_cl_metrics(met, sys, simParams, slsParams, slsOuts, s_a);
met
%% visualize
visualize_matrices(sys, slsOuts, slsOuts_alt, 20, 2);

%% calculate LQR cost
% note: have to separately save slsOutsCent
% also, we calculate the closed loop map up to T; have verified that
% elements after that are about 0

% slsParams_alt       = copy(slsParams);
% slsParams_alt.tFIR_ = Tc;
% [G, H] = calc_cl_map(sys, slsParams_alt, slsOuts_alt, simParams, slsParams_alt.tFIR_);
% calc_lqr_costs(slsOutsCent, G, H)

%% parameter sweep (Tc)
Tcs    = 2:slsParams.tFIR_+5;
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

%% check comm delay constraint violations (optional)
% useful to check when we have encouraged delay tolerance

delay = 1; % can adjust

[RZeros, MZeros] = get_delay_constraints(sys, Tc, delay);
Rc = slsOuts_alt.R_;
Mc = slsOuts_alt.M_;

RZeroViolations = 0;
MZeroViolations = 0;
totRZeroConstr  = 0;
totMZeroConstr  = 0;

for t=1:Tc
    RZeroViolations = RZeroViolations + sum(abs(Rc{t}(RZeros{t})) > tol);
    totRZeroConstr  = totRZeroConstr + sum(vec(RZeros{t} > 0)); 
    MZeroViolations = MZeroViolations + sum(abs(Mc{t}(MZeros{t})) > tol);
    totMZeroConstr  = totMZeroConstr + sum(vec(MZeros{t} > 0)); 
end

fprintf('%d of %d comm violations in R', RZeroViolations, totRZeroConstr);
fprintf('%d of %d comm violations in M', MZeroViolations, totMZeroConstr);
