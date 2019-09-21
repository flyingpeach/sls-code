%% choose the system you want to work with
setup1;

% for rank/zero conditions, try to match the precision of cvx_precision low
% http://cvxr.com/cvx/doc/solver.html#solver-precision
eps = 2.22e-16;
tol = eps.^(3/8);
close all;

settings = AltImplSettings(tol);

%% sandbox
% modes: ImplicitOpt, ExplicitOpt, Analytic, ApproxDrop, ApproxLeaky
settings.mode_      = AltImplMode.NullsOpt;
%settings.locality_  = 'encouraged';

settings.locality_ = 'strict';
settings.delay_ = 1;

Tc          = slsParams.tFIR_;

[RZeros, MZeros] = delay_constraints(sys, Tc, settings.delay_);
%settings.clDiffPen_ = 1e4;
%settings.relaxPct_  = 0.6;

slsOuts_alt = find_alt_impl(sys, slsParams, slsOuts, Tc, settings);
%%
s_a{1}  = slsOuts_alt;

met     = AltImplMetrics(tol, Tc);
met     = calc_mtx_metrics(met, sys, slsParams, slsOuts, s_a);
met     = calc_cl_metrics(met, sys, simParams, slsParams, slsOuts, s_a);
met

%% LQR/H2 costs

% inf-horizon centralized
[P,L,K] = dare(full(sys.A), full(sys.B2), eye(sys.Nx), eye(sys.Nu));

infH2Cost = 0;
for i=1:sys.Nx % H2 cost is sum of all costs from init condn
    x0    = zeros(sys.Nx, 1);
    x0(i) = 1;
    infH2Cost = infH2Cost + x0'*P*x0;
end

% centralized SLS (from setup1)
centH2Cost = 0;
for t=1:slsParams.tFIR_  
    centH2Cost = centH2Cost + norm(full([slsOuts.R_{t};slsOuts.M_{t}]), 'fro').^2;
end

% localized SLS
% for delay=1, 
locSLSH2Cost = 16.7383;

% centralized re-implementation
reImplH2Cost = 0;
tTotal = slsParams.tFIR_;

slsParams_alt = copy(slsParams);
slsParams_alt.tFIR_ = Tc;

[G, H] = calc_cl_map(sys, slsParams_alt, slsOuts_alt, simParams, tTotal);
for t=1:tTotal
    reImplH2Cost = reImplH2Cost + norm(full([G{t}; H{t}]), 'fro').^2;
end

% output all
[infH2Cost, centH2Cost, locSLSH2Cost, reImplH2Cost]

%% Internal stability
Tcs = 2:20;
specRadii = zeros(length(Tcs), 1);

for i=1:length(Tcs)
    specRadii(i) = check_int_stability(sys, Tcs(i), slsOuts_alts{i});
end

figure; hold on;
plot(Tcs, specRadii);
plot(Tcs, ones(length(Tcs),1), 'r');
xlabel('T_c'); ylabel('Spectral radius');

%%
visualize_matrices(sys, slsOuts, slsOuts_alt, Tc, 'all');

% check solver/feasibility statuses
disp(['Statuses:', print_statuses(sys, slsParams, slsOuts, s_a, tol)]);

%% check locality violation
% not used but keep for now
% slsParams_alt = copy(slsParams);
% slsParams_alt.tFIR_ = Tc;
% [RSupp, MSupp, ~] = get_localized_supports(sys, slsParams_alt);
% 
% Rc = slsOuts_alt.R_;
% Mc = slsOuts_alt.M_;
% 
% RLocViolations = 0;
% MLocViolations = 0;
% totalRSBZ = 0;
% totalMSBZ = 0;
% 
% for t=1:Tc
%     RShouldBZero   = Rc{t}(not(RSupp{t}));
%     RLocViolations = RLocViolations + sum(abs(RShouldBZero) > tol);
%     totalRSBZ      = totalRSBZ + sum(vec(not(RSupp{t}) > 0));
%     
%     MShouldBZero   = Mc{t}(not(MSupp{t}));
%     MLocViolations = MLocViolations + sum(abs(MShouldBZero) > tol);
%     totalMSBZ      = totalMSBZ + sum(vec(not(MSupp{t}) > 0));
% end
% 
% RLocViolations
% MLocViolations
% totalRSBZ
% totalMSBZ
%% find new impl over different Tcs
settings.mode_      = AltImplMode.ApproxLeaky;
settings.clDiffPen_ = 1e3;

Tcs    = 2:slsParams.tFIR_;
numTcs = length(Tcs);
slsOuts_alts = cell(numTcs, 1);
for idx=1:numTcs
    Tc = Tcs(idx);
    slsOuts_alts{idx} = find_alt_impl(sys, slsParams, slsOuts, Tc, settings);
end

scanH3 = slsOuts_alts;
%% find new impl over different approximations (ApproxDrop) 
Tc = round(slsParams.tFIR_/2);

settings       = AltImplSettings;
settings.mode_ = AltImplMode.ApproxDrop;
relaxPcts      = 0.05:0.05:1; 
numRelaxPcts   = length(relaxPcts);

slsOuts_alts   = cell(numRelaxPcts, 1);
for idx=1:numRelaxPcts
    settings.relaxPct_ = relaxPcts(idx);
    slsOuts_alts{idx}  = find_alt_impl(sys, slsParams, slsOuts, Tc, settings);
end

%% find new impl over different approximations (ApproxLeaky)
Tc = round(slsParams.tFIR_/4);

settings       = AltImplSettings;
settings.mode_ = AltImplMode.ApproxLeaky;
clDiffPens     = 1:6; % powers
numClDiffs     = length(clDiffPens);

slsOuts_alts   = cell(numClDiffs, 1);
for idx=1:numClDiffs
    settings.clDiffPen_ = 10^(clDiffPens(idx));
    slsOuts_alts{idx}  = find_alt_impl(sys, slsParams, slsOuts, Tc, settings);
end

%% check feasibility / solver statuses
disp(['Statuses:', print_statuses(sys, slsParams, slsOuts, slsOuts_alts, tol)]); 

%% plot stuff
% we might not want to plot all of slsOuts_alts, so this is the sandbox to
% adjust which slsOuts_alts to plot

sweepParamName = 'Tc';
%sweepParamName = 'relaxPct';
%sweepParamName = 'clDiffPen';

if strcmp(sweepParamName, 'Tc')
    xSeries = Tcs;
    xSize   = numTcs;
elseif strcmp(sweepParamName, 'relaxPct')
    xSeries = relaxPcts;
    xSize   = numRelaxPcts;
else
    xSeries = clDiffPens;
    xSize   = numClDiffs;    
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
    met = AltImplMetrics(tol, xPlot);
else
    met = AltImplMetrics(tol, Tc, sweepParamName, xPlot);
end

% calculate matrix-specific metrics
met = calc_mtx_metrics(met, sys, slsParams, slsOuts, slsOuts_altsPlot);

% calculate closed-loop metrics and plot select heat maps
met = calc_cl_metrics(met, sys, simParams, slsParams, slsOuts, slsOuts_altsPlot);

% plot metrics and save to file
savepath = 'C:\Users\Lisa\Desktop\caltech\research\implspace\tmp\';
plot_metrics(met, savepath);

% print again
disp(['Statuses:', print_statuses(sys, slsParams, slsOuts, slsOuts_alts, tol)]);  