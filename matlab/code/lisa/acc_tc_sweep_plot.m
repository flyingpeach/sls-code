%% set up system
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

for idx=1:numTcs
    solnSpaceSize = check_soln_space_size(sys, slsOutsCent, Tcs(idx), tol);
    disp([Tcs(idx), solnSpaceSize])
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
cStats = CtrllerStats(eps_nullsp, Tcs);
cStats = get_ctrller_stats(cStats, sys, slsOutsCent, csSweep);

fig1h  = figure(1);

subplot(1,3,1);
hold on;
p1=plot(Tcs, cStats.GcDiffs, '.-', 'MarkerSize', 16);
p2=plot(Tcs, cStats.HcDiffs, '.-', 'MarkerSize', 16);
title('Closed-loop differences'); legend('R', 'M');
xlabel('Tc'); ylabel('% difference'); 
xlim([1.5 25.5]); ylim([0 0.125]);
set(p1,'Color', [0 0.8 0]); set(p1,'LineWidth', 1.5);
set(p2,'Color', [0.8 0 0.8]); set(p2,'LineWidth', 1.5);

subplot(1,3,2);
hold on;
p3=plot(Tcs, cStats.IntSpecRadii_c, '.-', 'MarkerSize', 16);
p4=plot(Tcs, cStats.IntSpecRadiusOrig * ones(length(Tcs),1));
title('Spectral radii'); legend('new', 'original')
xlabel('Tc'); ylabel('spectral radius');
xlim([1.5 25.5]); ylim([0.1 1.2]);
set(p3,'Color', [0.8 0 0]); set(p3,'LineWidth', 1.5);
set(p4,'Color', [0 0 0.8]); set(p4,'LineWidth', 1.5);

subplot(1,3,3);
hold on;
p5=plot(Tcs, cStats.L1Norms, '.-', 'MarkerSize', 16);
p6=plot(Tcs, cStats.L1NormOrig * ones(length(Tcs), 1));
title('L1 norms'); legend('new', 'original')
xlabel('Tc'); ylabel('L1 norm');
xlim([1.5 25.5]);

set(p5,'Color', [0.8 0 0]); set(p5,'LineWidth', 1.5);
set(p6,'Color', [0 0 0.8]); set(p6,'LineWidth', 1.5);

x = 500; y = 200; width = 900; height = 200;
set(fig1h, 'Position', [x y width height]);

savepath = 'C:\Users\flyin\Desktop\caltech\research\sls controller design';
fn       = [savepath, '\all'];
fig1h.PaperUnits    = 'centimeters';
fig1h.PaperPosition = [0 0 22 6];
print(fn,'-dpng','-r300');