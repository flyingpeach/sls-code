Tcs = 2:25;
setup1;

approx_Tcs = 2:18;
approx_idx = 1:17;

exact_Tcs = 19:25;
exact_idx = 18:24;

% slsOuts_alts = slsOuts_alts_LS(approx_idx);
% slsOuts_alts(exact_idx) = slsOuts_alts_exact(exact_idx);
 
savepath = 'C:\Users\Lisa\Desktop\caltech\research\implspace\tmp\';

% not used
eps = 2.22e-16;
tol = eps.^(3/8);
zThresh = tol;

met = AltImplMetrics(zThresh, Tcs);
% calculate matrix-specific metrics
met = calc_mtx_metrics(met, sys, slsParams, slsOuts, slsOuts_alts);
% calculate closed-loop metrics and plot select heat maps
met = calc_cl_metrics(met, sys, simParams, slsParams, slsOuts, slsOuts_alts);

%%
% differences between CL maps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1h = figure(1);

subplot(1,3,1);
hold on;
p1=plot(Tcs, met.GcDiffs, '.-', 'MarkerSize', 16);
p2=plot(Tcs, met.HcDiffs, '.-', 'MarkerSize', 16);
title('Closed-loop differences');
legend('R', 'M');
ylabel('% difference'); xlabel('Tc');
xlim([1.5 25.5]);
ylim([0 0.125])

set(p1,'Color', [0 0.8 0]); set(p1,'LineWidth', 1.5);
set(p2,'Color', [0.8 0 0.8]); set(p2,'LineWidth', 1.5);

subplot(1,3,2);
hold on;
p3=plot(Tcs, met.IntSpecRadii_c, '.-', 'MarkerSize', 16);
p4=plot(Tcs, met.IntSpecRadiusOrig * ones(length(Tcs),1));
title('Spectral radii');
xlabel('Tc'); ylabel('spectral radius');
legend('new', 'original')
xlim([1.5 25.5]);
ylim([0.1 1.2]);

set(p3,'Color', [0.8 0 0]); set(p3,'LineWidth', 1.5);
set(p4,'Color', [0 0 0.8]); set(p4,'LineWidth', 1.5);

subplot(1,3,3);
hold on;
p5=plot(Tcs, met.L1Norms, '.-', 'MarkerSize', 16);
p6=plot(Tcs, met.L1NormOrig * ones(length(Tcs), 1));
title('L1 norms');
xlabel('Tc'); ylabel('L1 norm');
legend('new', 'original')
xlim([1.5 25.5]);

set(p5,'Color', [0.8 0 0]); set(p5,'LineWidth', 1.5);
set(p6,'Color', [0 0 0.8]); set(p6,'LineWidth', 1.5);

x = 500; y = 200;
width = 900; height = 200;
set(fig1h, 'Position', [x y width height]);

%%
fn = [savepath, '\all'];
fig1h.PaperUnits = 'centimeters';
fig1h.PaperPosition = [0 0 22 6];
print(fn,'-dpng','-r300');

% LQRCosts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Omitted for now; it actually is ... better?!
% figure; hold on;
% plot(xSeries, met.LQRCosts, 'o-');
% plot(xSeries, met.LQRCostOrig * ones(length(xSeries), 1));
% title('H2-LQR Costs');
% ylabel('H2-LQR Cost');
% legend('New LQR Costs', 'Original LQR Cost');
% savefig([savepath, '\lqr_costs.fig']);

%%
% check internal stability %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig2h = figure(2);
hold on;
p3=plot(Tcs, met.IntSpecRadii_c, '.-', 'MarkerSize', 16);
p4=plot(Tcs, met.IntSpecRadiusOrig * ones(length(Tcs),1));
title('Spectral radii of internal dynamics');
xlabel('Tc'); ylabel('spectral radius');
legend('new', 'original')
ylim([0.1 1.2]);

set(p3,'Color', [0.8 0 0]); set(p3,'LineWidth', 1.5);
set(p4,'Color', [0 0 0.8]); set(p4,'LineWidth', 1.5);

set(fig2h, 'Position', [x y width height]);

%%
fn = [savepath, '\spec_radii'];
fig2h.PaperUnits = 'centimeters';
fig2h.PaperPosition = [0 0 8 5];
print(fn,'-dpng','-r300') ;

%%
% compare L1 norms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig3h = figure(3); 
hold on;
p5=plot(Tcs, met.L1Norms, '.-', 'MarkerSize', 16);
p6=plot(Tcs, met.L1NormOrig * ones(length(Tcs), 1));
title('L1 norms');
xlabel('Tc'); ylabel('L1 norm');
legend('new', 'original')

set(p5,'Color', [0.8 0 0]); set(p5,'LineWidth', 1.5);
set(p6,'Color', [0 0 0.8]); set(p6,'LineWidth', 1.5);

set(fig3h, 'Position', [x y width height]);
%%
fn = [savepath, '\l1_norm'];
fig3h.PaperUnits = 'centimeters';
fig3h.PaperPosition = [0 0 8 5];
print(fn,'-dpng','-r300') ;

%% HEAT MAPS for example 2
simParams           = SimParams;
simParams.w_        = zeros(sys.Nx, 60); % disturbance
simParams.w_(1, 1) = 1;
simParams.openLoop_ = false;

simParams.tSim_     = 60;
[xorig1, uorig1] = simulate_system(sys, slsParams, slsOuts, simParams);
[xloc1, uloc1] = simulate_system(sys, slsParams, slsOuts_alt, simParams);

%simParams.tSim_     = 60;
[xvirt1, uvirt1] = simulate_system(sys, slsParams, slsOuts_virt, simParams);

simParams.w_(1, 1) = 0;
simParams.w_(5, 1) = 1;

%simParams.tSim_     = 40;
[xorig5, uorig5] = simulate_system(sys, slsParams, slsOuts, simParams);
[xloc5, uloc5] = simulate_system(sys, slsParams, slsOuts_alt, simParams);

%simParams.tSim_     = 60;
[xvirt5, uvirt5] = simulate_system(sys, slsParams, slsOuts_virt, simParams);

%% 
fig4h = figure(4);
logmin = -4; logmax = 0;

spaces = '                                                   ';
mySupTitle = ['Centralized' spaces 'Local' spaces 'Virtually Local'];
st = suptitle(mySupTitle);
set(st, 'FontSize', 11);
set(st, 'Position', [0.52 -0.02 0]);

subplot(2,6,1);
imagesc(log10(abs(xorig1))); caxis([logmin logmax]); 
title('log_1_0(|x|)'); ylabel('Space');
subplot(2,6,2);
imagesc(log10(abs(sys.B2*uorig1))); caxis([logmin logmax]);
title('log_1_0(|u|)');
set(gca, 'ytick',[]);

subplot(2,6,3);
imagesc(log10(abs(xloc1))); caxis([logmin logmax]); 
title('log_1_0(|x|)');
subplot(2,6,4);
imagesc(log10(abs(sys.B2*uloc1))); caxis([logmin logmax]);
title('log_1_0(|u|)');
set(gca, 'ytick',[])

subplot(2,6,5);
imagesc(log10(abs(xvirt1))); caxis([logmin logmax]); 
title('log_1_0(|x|)');
subplot(2,6,6);
imagesc(log10(abs(sys.B2*uvirt1))); caxis([logmin logmax]);
title('log_1_0(|u|)');
set(gca, 'ytick',[])

subplot(2,6,7);
imagesc(log10(abs(xorig5))); caxis([logmin logmax]); 
title('log_1_0(|x|)'); xlabel('Time'); ylabel('Space');
subplot(2,6,8);
imagesc(log10(abs(sys.B2*uorig5))); caxis([logmin logmax]);
title('log_1_0(|u|)'); xlabel('Time');
set(gca, 'ytick',[])

subplot(2,6,9);
imagesc(log10(abs(xloc5))); caxis([logmin logmax]); 
title('log_1_0(|x|)'); xlabel('Time');
subplot(2,6,10);
imagesc(log10(abs(sys.B2*uloc5))); caxis([logmin logmax]);
title('log_1_0(|u|)'); xlabel('Time');
set(gca, 'ytick',[])

subplot(2,6,11);
imagesc(log10(abs(xvirt5))); caxis([logmin logmax]); 
title('log_1_0(|x|)'); xlabel('Time');
subplot(2,6,12);
imagesc(log10(abs(sys.B2*uvirt5))); caxis([logmin logmax]);
title('log_1_0(|u|)'); xlabel('Time');
set(gca, 'ytick',[])

x = 100; y = 100;
width = 980; height = 345;
set(fig4h, 'Position', [x y width height]);

set(findall(gcf,'type','text'),'FontSize', 11)

h1=colorbar; colormap jet; 
set(h1, 'Position', [.92 .123 .02 .293])

h2=colorbar; colormap jet;
set(h2, 'Position', [.92 .596 .02 .293])

%%
fn = [savepath, '\heatmaps'];
fig4h.PaperUnits = 'centimeters';
fig4h.PaperPosition = [0 0 26 9.1];
print(fn,'-dpng','-r300') ;
