function plot_metrics(met, savepath)
% Plot trends of metrics from Rc, Mc
% Inputs
%     met       : AltImplMetrics containing calculated metrics
%     savepath  : where to save figures

if strcmp(met.sweepParamName, 'Tc')
    xSeries = met.Tcs;
else
    xSeries = met.sweepParams;
end
xLabel = met.sweepParamName;

% differences between CL maps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
subplot(2,1,1); hold on;
plot(xSeries, met.GcDiffs, 'o-');
plot(xSeries, met.GTcDiffs);
title('Normed diffs between CL maps from w to x');
ylabel('Normed difference');
legend('||Gc-R||_2', '||GTc-R||_2');

subplot(2,1,2); hold on;
plot(xSeries, met.HcDiffs, 'o-');
plot(xSeries, met.HTcDiffs);
title('Normed diffs between CL maps from w to u');
xlabel(xLabel); ylabel('Normed difference');
legend('||Hc-M||_2', '||HTc-M||_2');
savefig([savepath, '\cl_diff.fig']);

% compare L1 norms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; hold on;
plot(xSeries, met.L1Norms, 'o-');
plot(xSeries, met.L1NormsTc);
plot(xSeries, met.L1NormOrig * ones(length(xSeries), 1));
title('L1 norms of new implementation');
xlabel(xLabel); ylabel('L1-norm');
legend('L1 norms', 'L1 norms truncated', 'Original L1 norm');
savefig([savepath, '\l1_norms.fig']);

% compare number of nonzero entries %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
subplot(2,1,1); hold on;
plot(xSeries, met.RcNonzeros, 'o-');
plot(xSeries, met.RTcNonzeros);
plot(xSeries, met.RNonzero * ones(length(xSeries), 1));
title(sprintf('Entries of Rc/RTc/R > %0.1s', met.tol));
ylabel('# Nonzero entries');
legend('Rc', 'RTc', 'R');

subplot(2,1,2); hold on;
plot(xSeries, met.McNonzeros, 'o-');
plot(xSeries, met.MTcNonzeros);
plot(xSeries, met.MNonzero * ones(length(xSeries), 1));
title(sprintf('Entries of Mc/MTc/M > %0.1s', met.tol));
xlabel(xLabel); ylabel('# Nonzero entries');
legend('Mc', 'MTc', 'M');
savefig([savepath, '\nonzeros.fig']);
