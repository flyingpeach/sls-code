function met = calc_cl_metrics(met, sys, simParams, slsParams, slsOuts, slsOuts_alts, plotTcs, savepath)
% Calculate closed-loop metrics (i.e. compare x, u) of Rc, Mc and 
% corresponding RTc, MTc where RTc, MTc are R, M truncated to Tc elements
% also plots and saves heat maps for selected Tcs
%     met          : initialized AltImplMetrics that will be updated
% Inputs
%     sys          : LTISystem containing system matrices
%     simParams    : SimParams; parameters for the simulation
%     slsParams    : SLSParams containing parameters for original CL
%     slsOuts      : original closed loop system
%     slsOuts_alts : alternate CL implementations (one per different Tc)
%     plotTcs      : which Tc you want to plot a heatmap for
%     savepath     : where to save figures

numTcs = length(met.Tcs);

% original values
T = slsParams.tFIR_;
[x, u] = simulate_system(sys, slsParams, slsOuts, simParams);

% for simulation; will use Tc instead of T
slsParams_alt = copy(slsParams);

for i=1:numTcs
    Tc = met.Tcs(i);
    slsParams_alt.tFIR_ = Tc;

    % calculate CL metrics for Rc, Mc and plot heat maps %%%%%%%%%%%%%%%%%%
    [xNew, uNew]  = simulate_system(sys, slsParams_alt, slsOuts_alts{i}, simParams);
    met.xDiffs(i) = norm(xNew-x);
    met.uDiffs(i) = norm(uNew-u);

    if ismember(Tc, plotTcs)
        plot_heat_map(xNew, sys.B2*uNew, sprintf('New, Tc=%d', Tc));
        savefig([savepath, '\', sprintf('heatmap_Tc_%d.fig', Tc)]);
    end
    
    % calculate CL metrics for RTc, MTc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    if Tc < T
        [xTc, uTc] = simulate_system(sys, slsParams_alt, slsOuts, simParams);
        met.xTcDiffs(i) = norm(xTc-x);
        met.uTcDiffs(i) = norm(uTc-u);
    end
end

