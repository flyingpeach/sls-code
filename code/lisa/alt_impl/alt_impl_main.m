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

%% find new implementations Rc, Mc
Tcs     = [2:8];
numTcs  = length(Tcs);
plotTcs = []; % which Tcs you want to plot for

% slsOuts_alts{i} contains alternate implementation for Tc=Tcs(i)
slsOuts_alts = cell(numTcs, 1);

for i=1:numTcs
    Tc          = Tcs(i);
    slsOuts_alts{i} = find_alt_impl(sys, slsParams, slsOuts, Tc, 'approx');
end

%% calculate stats
zthresh = 1e-6; % below this value we'll count the value as a zero

% original values (renamed for easier reference)
T = slsParams.tFIR_;
R = slsOuts.R_; M = slsOuts.M_;

% stats for Rc, Mc
RDiffs     = zeros(numTcs, 1); MDiffs = zeros(numTcs, 1);
xDiffs     = zeros(numTcs, 1); uDiffs = zeros(numTcs, 1);
RcNonzeros = zeros(numTcs, 1); McNonzeros = zeros(numTcs, 1);
L1Norms    = zeros(numTcs, 1);

statusTxt = 'Statuses:'; % for writing out cvx statuses for each Tc

slsParams_alt = copy(slsParams); % for simulation; uses Tc instead of T

for i=1:numTcs
    Tc = Tcs(i); TMax = max(Tc, T);
    Rc = slsOuts_alts{i}.R_;
    Mc = slsOuts_alts{i}.M_;
   
    % make Rc, R equal length for easier comparisons
    for t=T+1:TMax % pad R, M with zeros if needed
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
        RcNonzeros(i) = RcNonzeros(i) + full(sum(vec(Rc{t} > zthresh)));
        McNonzeros(i) = McNonzeros(i) + full(sum(vec(Mc{t} > zthresh))); 
    end

    slsParams_alt.tFIR_ = Tc;
    [xNew, uNew] = simulate_system(sys, slsParams_alt, slsOuts_alts{i}, simParams);
    xDiffs(i)    = norm(xNew-xOld);
    uDiffs(i)    = norm(uNew-uOld);
    L1Norms(i)   = slsOuts_alts{i}.clnorm_;
    
    status    = slsOuts_alts{i}.solveStatus_;
    statusTxt = [statusTxt, char(10), sprintf('Tc=%d, %s', Tc, status)];

    if ismember(Tc, plotTcs)
        plot_heat_map(xNew, sys.B2*uNew, ['New, Tc=', int2str(Tc)]);
    end
end

% original stats
L1NormOrig = 0;
RNonzero   = 0; MNonzero  = 0; 

for t=1:T  
    L1NormOrig = L1NormOrig + norm([sys.C1, sys.D12]*[R{t}; M{t}], 1);
    RNonzero = RNonzero + full(sum(vec(R{t} > zthresh)));
    MNonzero = MNonzero + full(sum(vec(M{t} > zthresh)));
end

% reference stats
% our reference is truncated R, M (i.e. take only first Tc entries)
RTDiffs  = zeros(numTcs, 1); MTDiffs = zeros(numTcs, 1);
xTDiffs  = zeros(numTcs, 1); uTDiffs = zeros(numTcs, 1);
L1NormsT = zeros(numTcs, 1);

RTNonzeros = zeros(numTcs, 1);
MTNonzeros = zeros(numTcs, 1);

for i=1:numTcs
    Tc = Tcs(i);
    slsParams_alt.tFIR_ = Tc;
    
    for t=Tc+1:T
        RTDiffs(i) = RTDiffs(i) + norm(full(R{t}));
        MTDiffs(i) = MTDiffs(i) + norm(full(M{t}));
    end
   
   if Tc < T
       [xT, uT]   = simulate_system(sys, slsParams_alt, slsOuts, simParams);
   else
       xT = xOld; uT = uOld;
   end
   xTDiffs(i) = norm(xT-xOld);
   uTDiffs(i) = norm(uT-uOld);

   for t=1:Tc
       L1NormsT(i) = L1NormsT(i) + norm([sys.C1, sys.D12]*[R{t}; M{t}], 1);
       RTNonzeros(i) = RTNonzeros(i) + full(sum(vec(R{t} > zthresh)));
       MTNonzeros(i) = MTNonzeros(i) + full(sum(vec(M{t} > zthresh)));
   end
end

%% plot trends

savepath = 'C:\Users\Lisa\Desktop\caltech\research\implspace\tmp\';
disp(statusTxt);

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
savefig([savepath, 'mtxdiff.fig']);


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
savefig([savepath, 'cldiff.fig']);


figure; hold on; % compare L1 norms
plot(Tcs, L1Norms, 'o-');
plot(Tcs, L1NormsT, 'x-');
plot(Tcs, L1NormOrig * ones(numTcs, 1));
title('L1 norms of new implementation');
xlabel('Tc'); ylabel('L1-norm');
legend('L1 norms', 'L1 norms truncated', 'Original L1 norm');
savefig([savepath, 'l1normdiff.fig']);


figure; % compare # nonzero entries
subplot(2,1,1); hold on;
plot(Tcs, RcNonzeros, 'o-');
plot(Tcs, RTNonzeros, 'x-');
plot(Tcs, RNonzero * ones(numTcs, 1));
title(sprintf('Entries of Rc/RT/R > %0.1s', zthresh));
ylabel('# Nonzero entries');
legend('Rc', 'RT', 'R');

subplot(2,1,2); hold on;
plot(Tcs, McNonzeros, 'o-');
plot(Tcs, MTNonzeros, 'x-');
plot(Tcs, MNonzero * ones(numTcs, 1));
title(sprintf('Entries of Mc/MT/M > %0.1s', zthresh));
xlabel('Tc'); ylabel('# Nonzero entries');
legend('Mc', 'MT', 'M');
savefig([savepath, 'nonzeros.fig']);
