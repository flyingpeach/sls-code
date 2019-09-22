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
    met.L1NormOrig = met.L1NormOrig + norm([R{t}; M{t}], 1);
    met.RNonzero = met.RNonzero + sum(abs(vec(R{t})) > met.tol);
    met.MNonzero = met.MNonzero + sum(abs(vec(M{t})) > met.tol);
end
met.IntSpecRadiusOrig = check_int_stability(sys, T, R, M);

if strcmp(met.sweepParamName, 'Tc')
    numItems = numTcs;
else
    numItems = length(met.sweepParams);
    Tc       = met.Tcs;
end

for i=1:numItems
    if strcmp(met.sweepParamName, 'Tc')
        Tc = met.Tcs(i);
    end
    Rc = slsOuts_alts{i}.R_; Mc = slsOuts_alts{i}.M_;

    met.IntSpecRadii_c(i)  = check_int_stability(sys, Tc, Rc, Mc);
    met.IntSpecRadii_Tc(i) = check_int_stability(sys, Tc, R, M);
    
    for t=1:Tc
        % calculate metrics for Rc, Mc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        met.L1Norms(i) = met.L1Norms(i) + norm([Rc{t}; Mc{t}], 1);
        
        met.RcNonzeros(i) = met.RcNonzeros(i) + sum(abs(vec(Rc{t})) > met.tol);
        met.McNonzeros(i) = met.McNonzeros(i) + sum(abs(vec(Mc{t})) > met.tol);
        
        % calculate metrics for RTc, MTc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        met.L1NormsTc(i) = met.L1NormsTc(i) + norm([R{t}; M{t}], 1);
       
        met.RTcNonzeros(i) = met.RTcNonzeros(i) + sum(abs(vec(R{t})) > met.tol);
        met.MTcNonzeros(i) = met.MTcNonzeros(i) + sum(abs(vec(M{t})) > met.tol);
    end
end
