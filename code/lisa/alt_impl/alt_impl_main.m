%% choose the system you want to work with
setup1;

%% sandbox
Tc = 5;
slsOuts_alt = find_alt_impl(sys, slsParams, slsOuts, Tc);

slsParams_alt       = copy(slsParams);
slsParams_alt.tFIR_ = Tc;

[xNew, uNew] = simulate_system(sys, slsParams_alt, slsOuts_alt, simParams);

% Simulate but only use first Tc entries of original R, M
[xTrn, uTrn] = simulate_system(sys, slsParams_alt, slsOuts, simParams);

plot_heat_map(xOld, sys.B2*uOld, 'Original');
plot_heat_map(xTrn, sys.B2*uTrn, ['Truncated to Tc=', int2str(Tc)]);
plot_heat_map(xNew, sys.B2*uNew, ['New, Tc=', int2str(Tc)]);

L1NormOrig = 0;
L1NormTrunc = 0;
for t=1:slsParams.tFIR_    
    L1NormOrig = L1NormOrig + norm([sys.C1, sys.D12]*[slsOuts.R_{t}; slsOuts.M_{t}], 1);

    if t <= Tc % only count norms of first Tc terms
        L1NormTrunc = L1NormTrunc + norm([sys.C1, sys.D12]*[slsOuts.R_{t}; slsOuts.M_{t}], 1);
    end
end

L1NormOrig 
L1NormTrunc
L1NormNew = slsOuts_alt.clnorm_

%% find new implementations Rc, Mc; calculate stats
tFIR    = slsParams.tFIR_;
Tcs     = [2:tFIR, tFIR+1, tFIR+5];
plotTcs = [2, floor(tFIR/4), floor(tFIR/2), floor(tFIR*3/4)];

RDiffs  = zeros(length(Tcs), 1); MDiffs = zeros(length(Tcs), 1);
xDiffs  = zeros(length(Tcs), 1); uDiffs = zeros(length(Tcs), 1);
L1Norms = zeros(length(Tcs), 1);

slsParams_alt = copy(slsParams); % for simulation; uses Tc instead of tFIR

for i=1:length(Tcs)
    Tc = Tcs(i); TMax = max(Tc, tFIR);

    % find Rc, Mc (alternate implementations)
    slsOuts_alt = find_alt_impl(sys, slsParams, slsOuts, Tc);
    
    R  = slsOuts.R_; Rc = slsOuts_alt.R_;
    M  = slsOuts.M_; Mc = slsOuts_alt.M_;

    % calculate stats
    % make Rc, R equal length for easier comparisons
    for t=tFIR+1:TMax % pad R, M with zeros if needed
        R{t} = zeros(sys.Nx, sys.Nx);
        M{t} = zeros(sys.Nu, sys.Nx);
    end
    for t=Tc+1:TMax % pad Rc, Mc with zeros if needed
        Rc{t} = zeros(sys.Nx, sys.Nx);
        Mc{t} = zeros(sys.Nu, sys.Nx);
    end

    RDiff = 0; MDiff = 0;
    for t=1:TMax
        RDiff = RDiff + norm(full(R{t} - Rc{t}));
        MDiff = MDiff + norm(full(M{t} - Mc{t}));
    end
    RDiffs(i) = RDiff; MDiffs(i) = MDiff;

    slsParams_alt.tFIR_ = Tc;
    [xNew, uNew] = simulate_system(sys, slsParams_alt, slsOuts_alt, simParams);
    xDiffs(i)    = norm(xNew-xOld);
    uDiffs(i)    = norm(uNew-uOld);

    L1Norms(i) = slsOuts_alt.clnorm_;

    if ismember(Tc, plotTcs)
        plot_heat_map(xNew, sys.B2*uNew, ['New, Tc=', int2str(Tc)]);
    end
end

%% plot trends

figure; hold on;
plot(Tcs, RDiffs, 'o-');
plot(Tcs, MDiffs, 'o-');
title('Normed diffs between Rc/Mc, R/M');
xlabel('Tc'); ylabel('Normed difference');
legend('||Rc-R||_2', '||Mc-M||_2');

figure; hold on;
plot(Tcs, xDiffs, 'o-');
plot(Tcs, uDiffs, 'o-');
title('Normed diffs between CL responses');
xlabel('Tc'); ylabel('Normed difference');
legend('||xc-x||_2', '||uc-u||_2');

% L1 norm of original
L1NormOrig = 0;
for t=1:tFIR    
    L1NormOrig = L1NormOrig + norm([sys.C1, sys.D12]*[R{t}; M{t}], 1);
end

figure; hold on;
plot(Tcs, L1Norms, 'o-');
plot(Tcs, L1NormOrig * ones(length(Tcs), 1));
title('L1 norms of new implementation');
xlabel('Tc'); ylabel('L1-norm');
legend('L1 norms', 'Original L1 norm');