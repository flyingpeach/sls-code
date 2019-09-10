%% choose the system you want to work with
setup1;

%% sandbox
Tc = 2;
slsOuts_alt = find_alt_impl(sys, slsParams, slsOuts, Tc);

slsParams_alt       = copy(slsParams);
slsParams_alt.tFIR_ = Tc;

[xNew, uNew] = simulate_system(sys, slsParams_alt, slsOuts_alt, simParams);

plot_heat_map(xOld, sys.B2*uOld, 'Original');
plot_heat_map(xNew, sys.B2*uNew, ['New, Tc=', int2str(Tc)]);

%% find new implementations Rc, Mc; calculate stats
tFIR    = slsParams.tFIR_;
Tcs     = [2:25];
plotTcs = []; % which Tcs you want to plot for

thresh  = 1e-9; % below this value we'll count the value as a zero

RDiffs  = zeros(length(Tcs), 1); MDiffs = zeros(length(Tcs), 1);
xDiffs  = zeros(length(Tcs), 1); uDiffs = zeros(length(Tcs), 1);
L1Norms = zeros(length(Tcs), 1);

% count nonzero entries
RcNonzeros = zeros(length(Tcs), 1); 
McNonzeros = zeros(length(Tcs), 1);

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

    for t=1:TMax
        RDiffs(i) = RDiffs(i) + norm(full(R{t} - Rc{t}));
        MDiffs(i) = MDiffs(i) + norm(full(M{t} - Mc{t}));
        RcNonzeros(i) = RcNonzeros(i) + full(sum(vec(Rc{t} > thresh)));
        McNonzeros(i) = McNonzeros(i) + full(sum(vec(Mc{t} > thresh))); 
    end

    slsParams_alt.tFIR_ = Tc;
    [xNew, uNew] = simulate_system(sys, slsParams_alt, slsOuts_alt, simParams);
    xDiffs(i)    = norm(xNew-xOld);
    uDiffs(i)    = norm(uNew-uOld);

    L1Norms(i)  = slsOuts_alt.clnorm_;
    statuses{i} = slsOuts_alt.solveStatus_;
    
    if ismember(Tc, plotTcs)
        plot_heat_map(xNew, sys.B2*uNew, ['New, Tc=', int2str(Tc)]);
    end
end

%% calculate reference values

% original stats
L1NormOrig = 0;
RNonzero   = 0; MNonzero  = 0; 

for t=1:tFIR    
    L1NormOrig = L1NormOrig + norm([sys.C1, sys.D12]*[R{t}; M{t}], 1);
    RNonzero = RNonzero + full(sum(vec(R{t} > thresh)));
    MNonzero = MNonzero + full(sum(vec(M{t} > thresh)));
end

% reference stats
% our reference is truncated R, M (i.e. take only first Tc entries)
RTDiffs  = zeros(length(Tcs), 1); MTDiffs = zeros(length(Tcs), 1);
xTDiffs  = zeros(length(Tcs), 1); uTDiffs = zeros(length(Tcs), 1);
L1NormsT = zeros(length(Tcs), 1);

RTNonzeros = zeros(length(Tcs), 1);
MTNonzeros = zeros(length(Tcs), 1);

for i=1:length(Tcs)
    Tc = Tcs(i);
    slsParams_alt.tFIR_ = Tc;
    
    for t=Tc+1:tFIR
        RTDiffs(i) = RTDiffs(i) + norm(full(R{t}));
        MTDiffs(i) = MTDiffs(i) + norm(full(M{t}));
    end
   
   if Tc < tFIR
       [xT, uT]   = simulate_system(sys, slsParams_alt, slsOuts, simParams);
   else
       xT = xOld; uT = uOld;
   end
   xTDiffs(i) = norm(xT-xOld);
   uTDiffs(i) = norm(uT-uOld);

   for t=1:Tc
       L1NormsT(i) = L1NormsT(i) + norm([sys.C1, sys.D12]*[R{t}; M{t}], 1);
       RTNonzeros(i) = RTNonzeros(i) + full(sum(vec(R{t} > thresh)));
       MTNonzeros(i) = MTNonzeros(i) + full(sum(vec(M{t} > thresh)));
   end
end

%% plot trends

figure; % differences between the matrices
subplot(2,1,1); hold on;
plot(Tcs, RDiffs, 'o-');
plot(Tcs, RTDiffs, 'x-');
title('Normed diffs between R/Rc/RT');
ylabel('Normed difference');
legend('||Rc-R||_2', '||RT-R||_2');

subplot(2,1,2); hold on;
plot(Tcs, MDiffs, 'o-');
plot(Tcs, MTDiffs, 'x-');
title('Normed diffs between M/Mc/MT');
xlabel('Tc'); ylabel('Normed difference');
legend('||Mc-M||_2', '||MT-M||_2');


figure; % differences between CL state / inputs
subplot(2,1,1); hold on;
plot(Tcs, xDiffs, 'o-');
plot(Tcs, xTDiffs, 'x-');
title('Normed diffs between CL states');
ylabel('Normed difference');
legend('||xc-x||_2', '||xT-x||_2');

subplot(2,1,2); hold on;
plot(Tcs, uDiffs, 'o-');
plot(Tcs, uTDiffs, 'x-');
title('Normed diffs between CL inputs');
xlabel('Tc'); ylabel('Normed difference');
legend('||uc-u||_2', '||uT-u||_2');


figure; hold on; % compare L1 norms
plot(Tcs, L1Norms, 'o-');
plot(Tcs, L1NormsT, 'x-');
plot(Tcs, L1NormOrig * ones(length(Tcs), 1));
title('L1 norms of new implementation');
xlabel('Tc'); ylabel('L1-norm');
legend('L1 norms', 'L1 norms truncated', 'Original L1 norm');


figure; % compare # nonzero entries
subplot(2,1,1); hold on;
plot(Tcs, RcNonzeros, 'o-');
plot(Tcs, RTNonzeros, 'x-');
plot(Tcs, RNonzero * ones(length(Tcs), 1));
title(sprintf('Entries of Rc/RT/R > %0.1s', thresh));
ylabel('# Nonzero entries');
legend('Rc', 'RT', 'R');

subplot(2,1,2); hold on;
plot(Tcs, McNonzeros, 'o-');
plot(Tcs, MTNonzeros, 'x-');
plot(Tcs, MNonzero * ones(length(Tcs), 1));
title(sprintf('Entries of Mc/MT/M > %0.1s', thresh));
xlabel('Tc'); ylabel('# Nonzero entries');
legend('Mc', 'MT', 'M');
