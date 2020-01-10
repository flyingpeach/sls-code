function cStats = get_ctrller_stats(cStats, sys, slsOuts, ctrllers)
% Calculate stats for Rc, Mc and original R, M
%     cStats    : initialized CtrllerStats that will be updated
% Inputs
%     sys       : LTISystem containing system matrices
%     slsOuts   : SLSOutputs of original closed loop system
%     ctrllers  : array of Ctrller (one per value of swept parameter)

% original values
R = slsOuts.R_; M = slsOuts.M_;
T = length(R);

% calculate metrics for original R, M %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=1:T  
    cStats.L1NormOrig = cStats.L1NormOrig + norm([R{t}; M{t}], 1);
    cStats.RNonzero   = cStats.RNonzero + sum(abs(vec(R{t})) > cStats.tol);
    cStats.MNonzero   = cStats.MNonzero + sum(abs(vec(M{t})) > cStats.tol);
end
cStats.IntSpecRadiusOrig = check_int_stability(sys, R, M);

for t=1:T % TODO: this calculation (and other similar ones) assume that
          % output matrices are identity
    cStats.LQRCostOrig = cStats.LQRCostOrig + norm(full([R{t}; M{t}]), 'fro').^2;
end

% calculate metrics Rc, Mc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numItems = length(cStats.sweepParams);

% check CL map for extra time steps (expect zero)
extraT = 100; %TODO: a bit hacky. Include a sanity check?
tTotal = T + extraT;

% total CL norms (for normalization)
RTotalNorms = 0; MTotalNorms = 0;
for t=1:T
    RTotalNorms = RTotalNorms + norm(full(R{t}), 'fro');
    MTotalNorms = MTotalNorms + norm(full(M{t}), 'fro');
end

% zero pad R, M for easy comparison with CL
for t=T+1:tTotal
    R{t} = zeros(sys.Nx, sys.Nx);
    M{t} = zeros(sys.Nu, sys.Nx);
end

for i=1:numItems
    Rc = ctrllers{i}.Rc_; Mc = ctrllers{i}.Mc_; Tc = length(Rc);
    cStats.IntSpecRadii_c(i)  = check_int_stability(sys, Rc, Mc);

    for t=1:Tc
        cStats.L1Norms(i)    = cStats.L1Norms(i) + norm([Rc{t}; Mc{t}], 1);        
        cStats.RcNonzeros(i) = cStats.RcNonzeros(i) + sum(abs(vec(Rc{t})) > cStats.tol);
        cStats.McNonzeros(i) = cStats.McNonzeros(i) + sum(abs(vec(Mc{t})) > cStats.tol);        
    end
    
    % get CL map of Rc, Mc
    [Gc, Hc] = get_cl_map(sys, ctrllers{i}, tTotal);
    
    for t=1:tTotal
        cStats.LQRCosts(i) = cStats.LQRCosts(i) + norm(full([Gc{t}; Hc{t}]), 'fro').^2;
    end
    
    for t=1:tTotal % calculate differences
        cStats.GcDiffs(i) = cStats.GcDiffs(i) + norm(full(R{t} - Gc{t}), 'fro');
        cStats.HcDiffs(i) = cStats.HcDiffs(i) + norm(full(M{t} - Hc{t}), 'fro');
    end
    
    cStats.GcDiffs(i) = cStats.GcDiffs(i) / RTotalNorms;
    cStats.HcDiffs(i) = cStats.HcDiffs(i) / MTotalNorms; 
end