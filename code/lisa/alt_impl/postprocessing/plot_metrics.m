function plot_metrics(met, savepath)
% Plot trends of metrics from Rc, Mc
% Inputs
%     met       : AltImplMetrics containing calculated metrics
%     savepath  : where to save figures

% differences between the matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; 
subplot(2,1,1); hold on;
plot(met.Tcs, met.RDiffs, 'o-');
plot(met.Tcs, met.RTcDiffs, 'x-');
title('Normed diffs between R/Rc/RTc');
ylabel('Normed difference');
legend('||Rc-R||_2', '||RTc-R||_2');

subplot(2,1,2); hold on;
plot(met.Tcs, met.MDiffs, 'o-');
plot(met.Tcs, met.MTcDiffs, 'x-');
title('Normed diffs between M/Mc/MTc');
xlabel('Tc'); ylabel('Normed difference');
legend('||Mc-M||_2', '||MTc-M||_2');
savefig([savepath, '\mtx_diff.fig']);

% differences between CL states / inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
subplot(2,1,1); hold on;
plot(met.Tcs, met.xDiffs, 'o-');
plot(met.Tcs, met.xTcDiffs, 'x-');
title('Normed diffs between CL states');
ylabel('Normed difference');
legend('||xc-x||_2', '||xTc-x||_2');

subplot(2,1,2); hold on;
plot(met.Tcs, met.uDiffs, 'o-');
plot(met.Tcs, met.uTcDiffs, 'x-');
title('Normed diffs between CL inputs');
xlabel('Tc'); ylabel('Normed difference');
legend('||uc-u||_2', '||uTc-u||_2');
savefig([savepath, '\cl_diff.fig']);

% compare L1 norms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; hold on;
plot(met.Tcs, met.L1Norms, 'o-');
plot(met.Tcs, met.L1NormsTc, 'x-');
plot(met.Tcs, met.L1NormOrig * ones(length(met.Tcs), 1));
title('L1 norms of new implementation');
xlabel('Tc'); ylabel('L1-norm');
legend('L1 norms', 'L1 norms truncated', 'Original L1 norm');
savefig([savepath, '\l1_norms.fig']);

% compare number of nonzero entries %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
subplot(2,1,1); hold on;
plot(met.Tcs, met.RcNonzeros, 'o-');
plot(met.Tcs, met.RTcNonzeros, 'x-');
plot(met.Tcs, met.RNonzero * ones(length(met.Tcs), 1));
title(sprintf('Entries of Rc/RTc/R > %0.1s', met.zThresh));
ylabel('# Nonzero entries');
legend('Rc', 'RTc', 'R');

subplot(2,1,2); hold on;
plot(met.Tcs, met.McNonzeros, 'o-');
plot(met.Tcs, met.MTcNonzeros, 'x-');
plot(met.Tcs, met.MNonzero * ones(length(met.Tcs), 1));
title(sprintf('Entries of Mc/MTc/M > %0.1s', met.zThresh));
xlabel('Tc'); ylabel('# Nonzero entries');
legend('Mc', 'MTc', 'M');
savefig([savepath, '\nonzeros.fig']);
