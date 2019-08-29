clear; close all; clc; 

% specify system matrices
sys    = LTISystem;
sys.Nx = 50;

alpha = 0.4; rho = 1; actDens = 0.5;
generate_dbl_stoch_chain(sys, rho, actDens, alpha); % generate sys.A, sys.B2

sys.B1  = eye(sys.Nx); % used in simulation
sys.C1  = [speye(sys.Nx); sparse(sys.Nu, sys.Nx)]; % used in H2/HInf ctrl
sys.D12 = [sparse(sys.Nx, sys.Nu); speye(sys.Nu)];

% sls parameters
slsParams           = SLSParams;
slsParams.obj_      = Objective.H2;
slsParams.mode_     = SLSMode.ApproxDLocalized;
slsParams.tFIR_     = 10;
slsParams.actDelay_ = 1;
slsParams.d_        = 6;
slsParams.robCoeff_ = 1000;

% simulation parameters
simParams           = SimParams;
simParams.tSim_     = 25;
simParams.w_        = zeros(sys.Nx, simParams.tSim_); % disturbance
simParams.w_(floor(sys.Nx/2), 1) = 10;
simParams.openLoop_ = false;

cSpeeds = [2 1.5 1.4 1.3 1.2 1.1 1 0.9 0.8 0.7 0.6 0.5 0.4];
cSims   = [2 1 0.4]; % which comm speeds to simulate & plot

clnorms     = zeros(length(cSpeeds), 1); 
robustStabs = zeros(length(cSpeeds), 1);

for i=1:length(cSpeeds)
    slsParams.cSpeed_ = cSpeeds(i);
    slsOuts = state_fdbk_sls(sys, slsParams);

    if ismember(cSpeeds(i), cSims)
        [x, u] = simulate_system(sys, slsParams, slsOuts, simParams);
        plot_heat_map(x, sys.B2*u, ['Comms = ',num2str(cSpeeds(i))]);
    end

    clnorms(i)     = slsOuts.clnorm_;
    robustStabs(i) = slsOuts.robustStab_;    
end

figure;
p1=plot(cSpeeds, clnorms,'o-');
set(gca, 'xdir', 'reverse');
title([int2str(sys.Nx), ' Node Chain' ]);
xlabel('Comm Speed'); ylabel('Localized H_2-Norm Cost');

figure
p2=plot(cSpeeds,robustStabs,'o-');
set(gca, 'xdir', 'reverse');
title([int2str(sys.Nx), ' Node Chain' ]);
xlabel('Comm Speed'); ylabel('Stability Margin');

% These are the font settings from previous papers; uncomment as wanted
% set(gca,'FontSize',16,'fontWeight','bold')
% set(p1,'Color','red')
% set(p1,'LineWidth', 2)

