function plot_ctrller_stats(cstats, savepath)
% Plot trends of metrics from Rc, Mc
% Inputs
%     met       : AltImplMetrics containing calculated metrics
%     savepath  : where to save figures

xSeries = cstats.sweepParams;
xLabel  = cstats.sweepParamName;

% differences between CL maps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; hold on;
plot(xSeries, cstats.GcDiffs, 'o-');
plot(xSeries, cstats.HcDiffs, 'x-');
title('Normed diffs between CL maps');
xlabel(xLabel); ylabel('Normed difference');
legend('w to x', 'w to u');
savefig([savepath, '\cl_diff.fig']);

% LQR costs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; hold on;
plot(xSeries, cstats.LQRCosts, 'o-');
plot(xSeries, cstats.LQRCostOrig * ones(length(xSeries), 1));
title('H2-LQR Costs');
xlabel(xLabel); ylabel('H2-LQR Cost');
legend('New LQR Costs', 'Original LQR Cost');
savefig([savepath, '\lqr_costs.fig']);

% check internal stability %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; hold on;
plot(xSeries, cstats.IntSpecRadii_c, 'o-');
plot(xSeries, cstats.IntSpecRadiusOrig * ones(length(xSeries), 1));
title('Spectral radii of new implementation');
xlabel(xLabel); ylabel('Spectral radius');
legend('New', 'Original');
savefig([savepath, '\spec_radii.fig']);

% compare L1 norms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; hold on;
plot(xSeries, cstats.L1Norms, 'o-');
plot(xSeries, cstats.L1NormOrig * ones(length(xSeries), 1));
title('L1 norms of new implementation');
xlabel(xLabel); ylabel('L1-norm');
legend('L1 norms', 'Original L1 norm');
savefig([savepath, '\l1_norms.fig']);

% compare number of nonzero entries %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
subplot(2,1,1); hold on;
plot(xSeries, cstats.RcNonzeros, 'o-');
plot(xSeries, cstats.RNonzero * ones(length(xSeries), 1));
title(sprintf('Entries of Rc/R > %0.1s', cstats.tol));
ylabel('# Nonzero entries');
legend('Rc', 'R');

subplot(2,1,2); hold on;
plot(xSeries, cstats.McNonzeros, 'o-');
plot(xSeries, cstats.MNonzero * ones(length(xSeries), 1));
title(sprintf('Entries of Mc/M > %0.1s', cstats.tol));
xlabel(xLabel); ylabel('# Nonzero entries');
legend('Mc', 'M');
savefig([savepath, '\nonzeros.fig']);