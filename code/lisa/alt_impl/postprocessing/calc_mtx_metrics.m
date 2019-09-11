function met = calc_mtx_metrics(met, sys, slsParams, slsOuts, slsOuts_alts)
% Calculate metrics of Rc, Mc and corresponding RTc, MTc where RTc, MTc 
% are R, M truncated to Tc elements. These metrics are only related to
% the matrices themselves and are not specific to any simulation parameters
%     met          : initialized AltImplMetrics that will be updated
% Inputs
%     sys          : LTISystem containing system matrices
%     slsParams    : SLSParams containing parameters for original CL
%     slsOuts      : original closed loop system
%     slsOuts_alts : alternate CL implementations (one per different Tc)

numTcs = length(met.Tcs);

% original values
T = slsParams.tFIR_;
R = slsOuts.R_; M = slsOuts.M_;

% calculate metrics for original R, M %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=1:T  
    met.L1NormOrig = met.L1NormOrig + norm([sys.C1, sys.D12]*[R{t}; M{t}], 1);
    met.RNonzero = met.RNonzero + full(sum(vec(R{t} > met.zThresh)));
    met.MNonzero = met.MNonzero + full(sum(vec(M{t} > met.zThresh)));
end

for i=1:numTcs
    Tc = met.Tcs(i);
    Rc = slsOuts_alts{i}.R_; Mc = slsOuts_alts{i}.M_;
    
    % calculate metrics for Rc, Mc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    TMax = max(Tc, T);
   
    % make Rc, R equal length for easier comparisons
    for t=T+1:TMax % pad R, M with zeros if needed
        R{t} = zeros(sys.Nx, sys.Nx); M{t} = zeros(sys.Nu, sys.Nx);
    end
    for t=Tc+1:TMax % pad Rc, Mc with zeros if needed
        Rc{t} = zeros(sys.Nx, sys.Nx);Mc{t} = zeros(sys.Nu, sys.Nx);
    end

    for t=1:TMax
        met.L1Norms(i) = met.L1Norms(i) + norm([sys.C1, sys.D12]*[Rc{t}; Mc{t}], 1);
        
        met.RcNonzeros(i) = met.RcNonzeros(i) + full(sum(vec(Rc{t} > met.zThresh)));
        met.McNonzeros(i) = met.McNonzeros(i) + full(sum(vec(Mc{t} > met.zThresh)));

        met.RDiffs(i) = met.RDiffs(i) + norm(full(R{t} - Rc{t}));
        met.MDiffs(i) = met.MDiffs(i) + norm(full(M{t} - Mc{t}));
    end
    
    % calculate metrics for RTc, MTc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for t=1:Tc
       met.L1NormsTc(i) = met.L1NormsTc(i) + norm([sys.C1, sys.D12]*[R{t}; M{t}], 1);
       
       met.RTcNonzeros(i) = met.RTcNonzeros(i) + full(sum(vec(R{t} > met.zThresh)));
       met.MTcNonzeros(i) = met.MTcNonzeros(i) + full(sum(vec(M{t} > met.zThresh)));
    end

    for t=Tc+1:T
        met.RTcDiffs(i) = met.RTcDiffs(i) + norm(full(R{t}));
        met.MTcDiffs(i) = met.MTcDiffs(i) + norm(full(M{t}));
    end
end
