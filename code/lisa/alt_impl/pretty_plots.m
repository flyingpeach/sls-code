% Tcs = 2:25;
% 
% approx_Tcs = 2:18;
% approx_idx = 1:17;
% 
% exact_Tcs = 19:25;
% exact_idx = 18:24;
% 
% slsOuts_alts = slsOuts_alts_LS(approx_idx);
% slsOuts_alts(exact_idx) = slsOuts_alts_exact(exact_idx);

zThresh = tol;

met = AltImplMetrics(zThresh, Tcs);
% calculate matrix-specific metrics
met = calc_mtx_metrics(met, sys, slsParams, slsOuts, slsOuts_alts);
% calculate closed-loop metrics and plot select heat maps
met = calc_cl_metrics(met, sys, simParams, slsParams, slsOuts, slsOuts_alts);


xLabel = 'T_c';

% differences between CL maps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
subplot(2,1,1); hold on;
plot(approx_Tcs, met.GcDiffs(approx_idx), 'o-');
plot(exact_Tcs, met.GcDiffs(exact_idx), '*-');
title('Normed diffs between CL maps from w to x');
ylabel('Normed difference');
%legend('||Gc-R||_2', '||GTc-R||_2');

subplot(2,1,2); hold on;
plot(approx_Tcs, met.HcDiffs(approx_idx), 'o-');
plot(exact_Tcs, met.HcDiffs(exact_idx), '*-');
title('Normed diffs between CL maps from w to u');
xlabel(xLabel); ylabel('Normed difference');
%legend('||Hc-M||_2', '||HTc-M||_2');
savefig([savepath, '\cl_diff.fig']);

% LQRCosts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Omitted for now; it actually is ... better?!
% figure; hold on;
% plot(xSeries, met.LQRCosts, 'o-');
% plot(xSeries, met.LQRCostOrig * ones(length(xSeries), 1));
% title('H2-LQR Costs');
% ylabel('H2-LQR Cost');
% legend('New LQR Costs', 'Original LQR Cost');
% savefig([savepath, '\lqr_costs.fig']);

% check internal stability %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; hold on;
plot(approx_Tcs, met.IntSpecRadii_c(approx_idx), 'o-');
plot(exact_Tcs, met.IntSpecRadii_c(exact_idx), '*-');
plot(Tcs, met.IntSpecRadiusOrig * ones(length(Tcs),1));
plot(Tcs, ones(length(Tcs),1));
title('Spectral radii of new implementation');
xlabel(xLabel); ylabel('Spectral radius');
legend('Rc/Mc', 'Original');
savefig([savepath, '\spec_radii.fig']);

% compare L1 norms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; hold on;
plot(approx_Tcs, met.L1Norms(approx_idx), 'o-');
plot(exact_Tcs, met.L1Norms(exact_idx), '*-');
plot(Tcs, met.L1NormOrig * ones(length(Tcs), 1));
title('L1 norms of new implementation');
xlabel(xLabel); ylabel('L1-norm');
legend('L1 norms', 'Original L1 norm');
savefig([savepath, '\l1_norms.fig']);