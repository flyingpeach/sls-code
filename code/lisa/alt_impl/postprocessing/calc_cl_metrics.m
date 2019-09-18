function met = calc_cl_metrics(met, sys, simParams, slsParams, slsOuts, slsOuts_alts)
% Calculate closed-loop maps of Rc/Mc and compare to original R,M
% corresponding RTc, MTc where RTc, MTc are R, M truncated to Tc elements
%     met          : initialized AltImplMetrics that will be updated
% Inputs
%     sys          : LTISystem containing system matrices
%     simParams    : SimParams; parameters for the simulation
%     slsParams    : SLSParams containing parameters for original CL
%     slsOuts      : original closed loop system
%     slsOuts_alts : alternate CL implementations (one per different Tc)

numTcs = length(met.Tcs);

% original values
tFIR = slsParams.tFIR_;
R = slsOuts.R_; M = slsOuts.M_;

% check CL map a few steps after tFIR; expect zeros
extraT = 5;
tTotal = tFIR + extraT;

% zero pad R, M for easy comparisons
for t=tFIR+1:tTotal
    R{t} = zeros(sys.Nx, sys.Nx);
    M{t} = zeros(sys.Nu, sys.Nx);
end

% for simulation; will use Tc instead of T
slsParams_alt = copy(slsParams);
slsOutsTc     = copy(slsOuts);

if strcmp(met.sweepParamName, 'Tc')
    numItems = numTcs;
else
    numItems            = length(met.sweepParams);
    Tc                  = met.Tcs;
    slsParams_alt.tFIR_ = Tc;
end

for i=1:numItems
    if strcmp(met.sweepParamName, 'Tc')
        Tc                  = met.Tcs(i);
        slsParams_alt.tFIR_ = Tc;
    end

    % calculate CL map using Rc, Mc
    [Gc, Hc] = calc_cl_map(sys, slsParams_alt, slsOuts_alts{i}, simParams, tTotal);

    % calculate CL map using RTc, MTc
    % note: we do this even if Tc > tFIR to check numerical precision
    slsOutsTc.R_ = R(1:Tc);
    slsOutsTc.M_ = M(1:Tc);
    [GTc, HTc] = calc_cl_map(sys, slsParams_alt, slsOutsTc, simParams, tTotal);

    for t=1:tTotal % calculate differences
        met.GcDiffs(i) = met.GcDiffs(i) + norm(full(R{t} - Gc{t}));
        met.HcDiffs(i) = met.HcDiffs(i) + norm(full(M{t} - Hc{t}));

        met.GTcDiffs(i) = met.GTcDiffs(i) + norm(full(R{t} - GTc{t}));
        met.HTcDiffs(i) = met.HTcDiffs(i) + norm(full(M{t} - HTc{t}));
    end
end


% local functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [G, H] = calc_cl_map(sys, slsParams_alt, slsOuts_alt, simParams, tTotal)
simParams.tSim_ = tTotal + 1; % reactions to disturbance will start @ t=2

for t=1:tTotal
    G{t} = zeros(sys.Nx, sys.Nx); % CL map from w to x
    H{t} = zeros(sys.Nu, sys.Nx); % CL map from w to u
end

for i=1:sys.Nx
    simParams.w_ = zeros(sys.Nx, simParams.tSim_);
    simParams.w_(i, 1) = 1;
    [x, u]  = simulate_system(sys, slsParams_alt, slsOuts_alt, simParams);

    for t=1:tTotal
        G{t}(:,i) = x(:,t+1);
        H{t}(:,i) = u(:,t+1);
    end
end

