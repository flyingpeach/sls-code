%% generate stats for table
clear; close all; clc; 

% set up system / parameters
sys    = LTISystem();
sys.Nx = 10;

alpha = 0.4; rho = 1; actDens = 0.3;
generate_dbl_stoch_chain(sys, rho, actDens, alpha); % generate sys.A, sys.B2

sys.B1  = eye(sys.Nx);
sys.C1  = [speye(sys.Nx); sparse(sys.Nu, sys.Nx)];
sys.D12 = [sparse(sys.Nx, sys.Nu); speye(sys.Nu)];

slsParams       = SLSParams();
slsParams.T_    = 20;
slsParams.obj_  = Objective.H2;

% centralized SLS
slsParams.mode_ = SLSMode.Basic;
slsOutsCent     = state_fdbk_sls(sys, slsParams);

% infinite horizon
[P,L,K] = dare(full(sys.A), full(sys.B2), eye(sys.Nx), eye(sys.Nu));

infH2Cost = 0;
for i=1:sys.Nx % H2 cost is sum of all costs from init condition
    x0    = zeros(sys.Nx, 1);
    x0(i) = 1;
    infH2Cost = infH2Cost + x0'*P*x0;
end

% non-approx distributed SLS 
slsParams.mode_     = SLSMode.Delayed;
slsParams.actDelay_ = 1;
slsParams.cSpeed_   = 0.5;
slsParams.d_        = 3;
slsOutsDist         = state_fdbk_sls(sys, slsParams);

% approx distributed SLS
slsParams.approx_   = true;
slsParams.robCoeff_ = 1e3;
slsOutsVirtDist     = state_fdbk_sls(sys, slsParams);

% our method for distributed SLS
eps_base   = 2.22e-16;
eps_nullsp = eps_base.^(3/8);

cParams             = CtrllerParams();
cParams.mode_       = slsParams.mode_;

cParams.eps_nullsp_ = eps_nullsp;
cParams.T_          = slsParams.T_;

cParams.actDelay_   = slsParams.actDelay_;
cParams.cSpeed_     = slsParams.cSpeed_;
cParams.d_          = slsParams.d_;
cParams.CLDiffPen_  = 1e2;

ctrllerOurs = find_ctrller(sys, slsOutsCent, cParams);
cStatsOurs  = CtrllerStats(eps_nullsp);
csOurs{1}   = ctrllerOurs;
cStatsOurs  = get_ctrller_stats(cStatsOurs, sys, slsOutsCent, csOurs);

% Make a "Controller" out of slsOutsDist
csDist         = Ctrller();
csDist.Rc_     = slsOutsDist.R_;
csDist.Mc_     = slsOutsDist.M_;
csDists{1}     = csDist;
cStatsDist     = CtrllerStats(eps_nullsp);

% Note: might have to comment out this line if the distributed
% problem was so infeasible that there's NaNs in the matrix
cStatsDist     = get_ctrller_stats(cStatsDist, sys, slsOutsCent, csDists);

% Make a "Controller" out of slsOutsVirtDist
csVirtDist     = Ctrller();
csVirtDist.Rc_ = slsOutsVirtDist.R_; 
csVirtDist.Mc_ = slsOutsVirtDist.M_;
csVirtDists{1} = csVirtDist;
cStatsVirtDist = CtrllerStats(eps_nullsp);
cStatsVirtDist = get_ctrller_stats(cStatsVirtDist, sys, slsOutsCent, csVirtDists);

disp(sprintf('Centralized  : cost=%.3f, rad=%.3f, l1norm=%.3f', ...
             cStatsOurs.LQRCostOrig / infH2Cost, ...
             cStatsOurs.IntSpecRadiusOrig, ...
             cStatsOurs.L1NormOrig))

disp(sprintf('Distr SLS    : cost=%.3f, rad=%.3f, l1norm=%.3f', ...
             cStatsDist.LQRCosts / infH2Cost, ...
             cStatsDist.IntSpecRadii_c, ...
             cStatsDist.L1Norms))         
         
disp(sprintf('VirtDist SLS : cost=%.3f, rad=%.3f, l1norm=%.3f', ...
             cStatsVirtDist.LQRCosts / infH2Cost, ...
             cStatsVirtDist.IntSpecRadii_c, ...
             cStatsVirtDist.L1Norms))

disp(sprintf('Our method   : cost=%.3f, rad=%.3f, l1norm=%.3f', ...
             cStatsOurs.LQRCosts / infH2Cost, ...
             cStatsOurs.IntSpecRadii_c, ...
             cStatsOurs.L1Norms))

%% plot closed-loop behaviour
simParams           = SimParams();
simParams.w_        = zeros(sys.Nx, 60);
simParams.openLoop_ = false;
simParams.tSim_     = 60;

% simulate with disturbance @ node 1
simParams.w_(1, 1) = 1;

[xcent1, ucent1] = simulate_system(sys, simParams, slsOutsCent.R_, slsOutsCent.M_);
[xvirt1, uvirt1] = simulate_system(sys, simParams, slsOutsVirtDist.R_, slsOutsVirtDist.M_);
[xours1, uours1] = simulate_system(sys, simParams, ctrllerOurs.Rc_, ctrllerOurs.Mc_);

% simulate with disturbance @ node 5
simParams.w_(1, 1) = 0;
simParams.w_(5, 1) = 1;

[xcent5, ucent5] = simulate_system(sys, simParams, slsOutsCent.R_, slsOutsCent.M_);
[xvirt5, uvirt5] = simulate_system(sys, simParams, slsOutsVirtDist.R_, slsOutsVirtDist.M_);
[xours5, uours5] = simulate_system(sys, simParams, ctrllerOurs.Rc_, ctrllerOurs.Mc_);

fig4h = figure(4);
logmin = -4; logmax = 0;

spaces1 = '                                      ';
spaces2 = '                                          ';
mySupTitle = ['Centralized' spaces1 'Virtually Distr' spaces2 'Distr'];
st = suptitle(mySupTitle);
set(st, 'FontSize', 13);
set(st, 'Position', [0.52 -0.02 0]);

subplot(2,6,1);
imagesc(log10(abs(xcent1))); caxis([logmin logmax]); 
title('log_1_0(|x|)'); ylabel('Space');
subplot(2,6,2);
imagesc(log10(abs(sys.B2*ucent1))); caxis([logmin logmax]);
title('log_1_0(|u|)');
set(gca, 'ytick',[]);

subplot(2,6,3);
imagesc(log10(abs(xvirt1))); caxis([logmin logmax]); 
title('log_1_0(|x|)');
subplot(2,6,4);
imagesc(log10(abs(sys.B2*uvirt1))); caxis([logmin logmax]);
title('log_1_0(|u|)');
set(gca, 'ytick',[])

subplot(2,6,5);
imagesc(log10(abs(xours1))); caxis([logmin logmax]); 
title('log_1_0(|x|)');
subplot(2,6,6);
imagesc(log10(abs(sys.B2*uours1))); caxis([logmin logmax]);
title('log_1_0(|u|)');
set(gca, 'ytick',[])

subplot(2,6,7);
imagesc(log10(abs(xcent5))); caxis([logmin logmax]); 
title('log_1_0(|x|)'); xlabel('Time'); ylabel('Space');
subplot(2,6,8);
imagesc(log10(abs(sys.B2*ucent5))); caxis([logmin logmax]);
title('log_1_0(|u|)'); xlabel('Time');
set(gca, 'ytick',[])

subplot(2,6,9);
imagesc(log10(abs(xvirt5))); caxis([logmin logmax]); 
title('log_1_0(|x|)'); xlabel('Time');
subplot(2,6,10);
imagesc(log10(abs(sys.B2*uvirt5))); caxis([logmin logmax]);
title('log_1_0(|u|)'); xlabel('Time');
set(gca, 'ytick',[])

subplot(2,6,11);
imagesc(log10(abs(xours5))); caxis([logmin logmax]); 
title('log_1_0(|x|)'); xlabel('Time');
subplot(2,6,12);
imagesc(log10(abs(sys.B2*uours5))); caxis([logmin logmax]);
title('log_1_0(|u|)'); xlabel('Time');
set(gca, 'ytick',[])

x = 100; y = 100; width = 980; height = 345;
set(fig4h, 'Position', [x y width height]);
set(findall(gcf,'type','text'),'FontSize', 12)

h1=colorbar; colormap jet; 
set(h1, 'Position', [.92 .123 .02 .293])
h2=colorbar; colormap jet;
set(h2, 'Position', [.92 .596 .02 .293])

savepath = 'C:\Users\flyin\Desktop\caltech\research\sls controller design';
fn       = [savepath, '\heatmaps'];
fig4h.PaperUnits    = 'centimeters';
fig4h.PaperPosition = [0 0 26 9.1];
print(fn,'-dpng','-r300') ;
