%% set up system
clear; close all; clc; 

sys    = LTISystem();
sys.Nx = 10;

alpha = 0.4; rho = 1; actDens = 0.3;
generate_dbl_stoch_chain(sys, rho, actDens, alpha); % generate sys.A, sys.B2

sys.B1  = eye(sys.Nx);
sys.C1  = [speye(sys.Nx); sparse(sys.Nu, sys.Nx)];
sys.D12 = [sparse(sys.Nx, sys.Nu); speye(sys.Nu)];

% sls paracStatsers
slsParams       = SLSParams();
slsParams.T_    = 20;
slsParams.obj_  = Objective.H2;
slsParams.mode_ = SLSMode.Basic;

slsOutsCent = state_fdbk_sls(sys, slsParams);
eps_nullsp  = 1e-2;

%% check solutions space sizes
Tcs    = 2:25;
numTcs = length(Tcs);
tol    = eps_nullsp;

eps_base = 1e-7;
for idx=1:numTcs
    F          = get_ctrller_constraint(sys, slsOutsCent, Tcs(idx));
    F2         = F(:,sys.Nx+1:end);
    solnExist  = ( rank(F, eps_base) == rank(F2, eps_base) );
    nullspSize = size(get_nullsp(F2, tol), 2);

    disp([Tcs(idx), solnExist, nullspSize]);
end

%% find controllers
cParams             = CtrllerParams();
cParams.mode_       = SLSMode.Basic;
cParams.eps_nullsp_ = eps_nullsp;
cParams.T_          = slsParams.T_;

csSweep  = cell(numTcs, 1);
for idx=1:numTcs
    cParams.T_ = Tcs(idx);
    csSweep{idx}  = find_ctrller(sys, slsOutsCent, cParams);
end

%% plot sweep over Tcs
cStats = CtrllerStats(eps_nullsp, 'Tc', Tcs);
cStats = get_ctrller_stats(cStats, sys, slsOutsCent, csSweep);

fig1h  = figure(1);

subplot(1,4,1);
hold on;
p1=plot(Tcs, cStats.GcDiffs, '.-', 'MarkerSize', 16);
p2=plot(Tcs, cStats.HcDiffs, '.-', 'MarkerSize', 16);
title('Closed-loop differences'); legend('R', 'M');
xlabel('Tc'); ylabel('% difference'); 
xlim([1.5 25.5]);
set(p1,'Color', [0 0.8 0]); set(p1,'LineWidth', 1.5);
set(p2,'Color', [0.8 0 0.8]); set(p2,'LineWidth', 1.5);

subplot(1,4,2);
hold on;
p3=plot(Tcs, cStats.IntSpecRadii_c, '.-', 'MarkerSize', 16);
p4=plot(Tcs, cStats.IntSpecRadiusOrig * ones(length(Tcs),1));
title('Spectral radii'); legend('new', 'original')
xlabel('Tc'); ylabel('spectral radius');
xlim([1.5 25.5]); 
set(p3,'Color', [0.8 0 0]); set(p3,'LineWidth', 1.5);
set(p4,'Color', [0 0 0.8]); set(p4,'LineWidth', 1.5);

subplot(1,4,3);
hold on;
p5=plot(Tcs, cStats.L1Norms, '.-', 'MarkerSize', 16);
p6=plot(Tcs, cStats.L1NormOrig * ones(length(Tcs), 1));
title('L1 norms'); legend('new', 'original')
xlabel('Tc'); ylabel('L1 norm');
xlim([1.5 25.5]);
set(p5,'Color', [0.8 0 0]); set(p5,'LineWidth', 1.5);
set(p6,'Color', [0 0 0.8]); set(p6,'LineWidth', 1.5);

subplot(1,4,4);
hold on;
p7=plot(Tcs, cStats.LQRCosts, '.-', 'MarkerSize', 16);
p8=plot(Tcs, cStats.LQRCostOrig * ones(length(Tcs), 1));
title('LQR costs'); legend('new', 'original')
xlabel('Tc'); ylabel('LQR cost');
xlim([1.5 25.5]);
set(p7,'Color', [0.8 0 0]); set(p7,'LineWidth', 1.5);
set(p8,'Color', [0 0 0.8]); set(p8,'LineWidth', 1.5);

x = 0; y = 200; width = 1050; height = 200;
set(fig1h, 'Position', [x y width height]); % display only
savepath = 'C:\Users\flyin\Desktop\caltech\research\sls controller design';
fn       = [savepath, '\all'];
fig1h.PaperUnits    = 'centimeters';
fig1h.PaperPosition = [0 0 30 6]; % print properties
print(fn,'-dpng','-r300');