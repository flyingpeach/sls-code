function met = get_ctrller_stats(met, sys, simParams, slsParams, slsOuts, ctrllers)
% Calculate stats for Rc, Mc and original R, M
%     met       : initialized CtrllerStats that will be updated
% Inputs
%     sys       : LTISystem containing system matrices
%     simParams : SimParams; parameters for the simulation
%     slsParams : SLSParams containing parameters for original CL
%     slsOuts   : SLSOutputs of original closed loop system
%     ctrllers  : array of Ctrller (one per value of swept parameter)

numTcs = length(met.Tcs);

% original values
tFIR = slsParams.tFIR_;
R = slsOuts.R_; M = slsOuts.M_;

% calculate metrics for original R, M %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=1:tFIR  
    met.L1NormOrig = met.L1NormOrig + norm([R{t}; M{t}], 1);
    met.RNonzero   = met.RNonzero + sum(abs(vec(R{t})) > met.tol);
    met.MNonzero   = met.MNonzero + sum(abs(vec(M{t})) > met.tol);
end
met.IntSpecRadiusOrig = check_int_stability(sys, tFIR, R, M);

for t=1:tFIR
    met.LQRCostOrig = met.LQRCostOrig + norm(full([R{t}; M{t}]), 'fro').^2;
end

% calculate metrics Rc, Mc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
slsParams_alt = copy(slsParams); % for simulation; will use Tc instead of T

if strcmp(met.sweepParamName, 'Tc')
    numItems = numTcs;
else
    numItems            = length(met.sweepParams);
    Tc                  = met.Tcs;
    slsParams_alt.tFIR_ = Tc;
end

% check CL map for extra time steps (expect zero)
extraT = 100;
tTotal = tFIR + extraT;

% total CL norms (for normalization)
RTotalNorms = 0; MTotalNorms = 0;
for t=1:tFIR
    RTotalNorms = RTotalNorms + norm(full(R{t}));
    MTotalNorms = MTotalNorms + norm(full(M{t}));
end

% zero pad R, M for easy comparison with CL
for t=tFIR+1:tTotal
    R{t} = zeros(sys.Nx, sys.Nx);
    M{t} = zeros(sys.Nu, sys.Nx);
end

for i=1:numItems
    if strcmp(met.sweepParamName, 'Tc')
        Tc                  = met.Tcs(i);
        slsParams_alt.tFIR_ = Tc;
    end
    Rc = ctrllers{i}.Rc_; Mc = ctrllers{i}.Mc_;
    met.IntSpecRadii_c(i)  = check_int_stability(sys, Tc, Rc, Mc);

    for t=1:Tc
        met.L1Norms(i)    = met.L1Norms(i) + norm([Rc{t}; Mc{t}], 1);        
        met.RcNonzeros(i) = met.RcNonzeros(i) + sum(abs(vec(Rc{t})) > met.tol);
        met.McNonzeros(i) = met.McNonzeros(i) + sum(abs(vec(Mc{t})) > met.tol);        
    end
    
    % get CL map of Rc, Mc
    [Gc, Hc] = get_cl_map(sys, slsParams_alt, ctrllers{i}, simParams, tTotal);
    
    for t=1:tTotal
        met.LQRCosts(i) = met.LQRCosts(i) + norm(full([Gc{t}; Hc{t}]), 'fro').^2;
    end
    
    for t=1:tTotal % calculate differences
        met.GcDiffs(i) = met.GcDiffs(i) + norm(full(R{t} - Gc{t}));
        met.HcDiffs(i) = met.HcDiffs(i) + norm(full(M{t} - Hc{t}));
    end
    
    met.GcDiffs(i) = met.GcDiffs(i) / RTotalNorms;
    met.HcDiffs(i) = met.HcDiffs(i) / MTotalNorms; 
end