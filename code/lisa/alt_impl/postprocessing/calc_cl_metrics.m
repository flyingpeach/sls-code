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

% original LQR cost
for t=1:tFIR
    met.LQRCostOrig = met.LQRCostOrig + norm(full([R{t}; M{t}]), 'fro').^2;
end

% check CL map a few steps after tFIR; expect zeros
extraT = 20;
tTotal = tFIR + extraT;

% zero pad R, M for easy comparisons
RTotalNorms = 0;
MTotalNorms = 0;

for t=1:tFIR
    RTotalNorms = RTotalNorms + norm(full(R{t}));
    MTotalNorms = MTotalNorms + norm(full(M{t}));
end

for t=tFIR+1:tTotal
    R{t} = zeros(sys.Nx, sys.Nx);
    M{t} = zeros(sys.Nu, sys.Nx);
end

% for simulation; will use Tc instead of T
slsParams_alt = copy(slsParams);
%slsOutsTc     = copy(slsOuts);

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
    for t=1:tTotal
        met.LQRCosts(i) = met.LQRCosts(i) + norm(full([Gc{t}; Hc{t}]), 'fro').^2;
    end   
    
    % calculate CL map using RTc, MTc
    % note: we do this even if Tc > tFIR to check numerical precision
%     slsOutsTc.R_ = R(1:Tc);
%     slsOutsTc.M_ = M(1:Tc);
%     [GTc, HTc] = calc_cl_map(sys, slsParams_alt, slsOutsTc, simParams, tTotal);

    for t=1:tTotal % calculate differences
        met.GcDiffs(i) = met.GcDiffs(i) + norm(full(R{t} - Gc{t}));
        met.HcDiffs(i) = met.HcDiffs(i) + norm(full(M{t} - Hc{t}));

%         met.GTcDiffs(i) = met.GTcDiffs(i) + norm(full(R{t} - GTc{t}));
%         met.HTcDiffs(i) = met.HTcDiffs(i) + norm(full(M{t} - HTc{t}));
    end
    
    met.GcDiffs(i) = met.GcDiffs(i) / RTotalNorms;
    met.HcDiffs(i) = met.HcDiffs(i) / MTotalNorms;
%     met.GTcDiffs(i) = met.GTcDiffs(i) / RTotalNorms;
%     met.HTcDiffs(i) = met.HTcDiffs(i) / MTotalNorms;
    
    
end
