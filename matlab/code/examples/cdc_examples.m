% note: run sequentially; some params from later sections depend on
% params specified in previous sections
clear; close all; clc; 

% specify system matrices
sys    = LTISystem;
sys.Nx = 50;

alpha = 0.4; rho = 1; actDens = 1;
generate_dbl_stoch_chain(sys, rho, actDens, alpha); % generate sys.A, sys.B2

sys.B1  = eye(sys.Nx); % used in simulation
sys.C1  = [speye(sys.Nx); sparse(sys.Nu, sys.Nx)]; % used in H2/HInf ctrl
sys.D12 = [sparse(sys.Nx, sys.Nu); speye(sys.Nu)];

% sls parameters
slsParams           = SLSParams;
slsParams.obj_      = Objective.H2;
slsParams.tFIR_     = 10;

% simulation parameters
simParams           = SimParams;
simParams.w_        = zeros(sys.Nx, 100);
simParams.w_(floor(sys.Nx/2), 1) = 10;

%% open loop and centralized solutions
slsParams.mode_ = SLSMode.Basic;
slsOutsCent     = state_fdbk_sls(sys, slsParams);

% simulate open and closed loop
simParams.openLoop_ = true; 
simParams.tSim_ = 90;
[xOpen, uOpen]  = simulate_system(sys, slsParams, slsOutsCent, simParams);

simParams.openLoop_ = false; 
simParams.tSim_ = 25;
[xCent, uCent]  = simulate_system(sys, slsParams, slsOutsCent, simParams);

plot_heat_map(xOpen, sys.B2*uOpen, 'Open loop');
plot_heat_map(xCent, sys.B2*uCent, [int2str(sys.Nx), ' Node Centralized Solution']);

%% varying horizon lengths (tFIR)
tFIRs   = [3, 4, 5, 6, 7, 8];
tPrints = [4, 8];

slsParams.mode_     = SLSMode.DLocalized;
slsParams.actDelay_ = 1;
slsParams.cSpeed_   = 2;
slsParams.d_        = 6;
simParams.tSim_     = 16;
clnorms             = zeros(length(tFIRs), 1);
for i=1:length(tFIRs)
    slsParams.tFIR_ = tFIRs(i);
    slsOutsT        = state_fdbk_sls(sys, slsParams);
    clnorms(i)      = slsOutsT.clnorm_;
    
    if ismember(tFIRs(i), tPrints)   
        [xT, uT] = simulate_system(sys, slsParams, slsOutsT, simParams);
        plot_heat_map(xT, sys.B2*uT, ['tFIR = ', int2str(tFIRs(i))]);
    end
end

figure;
p1=plot(tFIRs, clnorms, 'o-');
set(gca, 'xdir', 'reverse');
title([int2str(sys.Nx), ' Node Chain']);
xlabel('FIR Horizon T'); ylabel('Localized H_2-Norm Cost');

% These are the font settings from previous papers; uncomment as wanted
% also applicable to other sections on this plot; just set the plot name
% (i.e. p1) appropriately

% set(gca,'FontSize',16,'fontWeight','bold');
% set(p1,'Color','red'); set(p1,'LineWidth', 2);

%% varying actuation densities (actDens)
actDenss  = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];
actPrints = [0.2, 0.3, 0.7, 1];

slsParams.tFIR_ = 15;    
simParams.tSim_ = 45;
clnorms         = zeros(length(actDenss), 1);
for i = 1:length(actDenss)
    generate_dbl_stoch_chain(sys, rho, actDenss(i), alpha); % update A, B2
    % Nu changed, so have to change these as well
    sys.C1  = [speye(sys.Nx); sparse(sys.Nu, sys.Nx)];
    sys.D12 = [sparse(sys.Nx, sys.Nu); speye(sys.Nu)];

    if actDenss(i) == 0.2
        % in this case density is too low; use approx instead of direct
        slsParams.mode_     = SLSMode.ApproxDLocalized;
        slsParams.robCoeff_ = 1000;
    else
        slsParams.mode_ = SLSMode.DLocalized;
    end

    slsOutsAct = state_fdbk_sls(sys, slsParams);
    clnorms(i) = slsOutsAct.clnorm_;
    
    if ismember(actDenss(i), actPrints)
        [xAct, uAct] = simulate_system(sys, slsParams, slsOutsAct, simParams);
        plot_heat_map(xAct, sys.B2*uAct, ['Actuation = ', num2str(actDenss(i))]);
    end
end
   
figure;
p2=plot(actDenss, clnorms, 'o-');
set(gca, 'xdir', 'reverse');
title([int2str(sys.Nx), ' Node Chain']);
xlabel('Actuation Density'); ylabel('Localized H_2-Norm Cost');

%% varying locality constraints (d)
ds      = [2, 3, 4, 5, 6, 8, 10];
dPrints = [2, 4, 10];

alpha = 0.4; rho = 1; actDens = 1;
generate_dbl_stoch_chain(sys, rho, actDens, alpha);
% Nu changed, so have to change these as well
sys.C1  = [speye(sys.Nx); sparse(sys.Nu, sys.Nx)];
sys.D12 = [sparse(sys.Nx, sys.Nu); speye(sys.Nu)];

slsParams.tFIR_ = 10;
simParams.tSim_ = 16;
clnorms         = zeros(length(ds), 1);
for i=1:length(ds)
    slsParams.d_ = ds(i);
    if ds(i) == 2
        % in this case locality is too strict; use approx instead of direct
        slsParams.mode_     = SLSMode.ApproxDLocalized;
        slsParams.robCoeff_ = 1000;
    else
        slsParams.mode_ = SLSMode.DLocalized;
    end

    slsOutsD   = state_fdbk_sls(sys, slsParams);
    clnorms(i) = slsOutsD.clnorm_;
    
    if ismember(ds(i), dPrints)
        [xD, uD] = simulate_system(sys, slsParams, slsOutsD, simParams);
        plot_heat_map(xD, sys.B2*uD, ['Locality = ', int2str(ds(i))]);
    end   
end

figure;
p3=plot(ds, clnorms, 'o-');
set(gca, 'xdir', 'reverse');
title([int2str(sys.Nx), ' Node Chain'])
xlabel('d-hops'); ylabel('Localized H_2-Norm Cost');

%% varying communication speeds (cSpeed)
cSpeeds = [1, 1.25, 1.5, 1.75, 2, 3, 4];
cPrints = [1, 1.75, 4];

slsParams.d_    = 8;
simParams.tSim_ = 25;
clnorms         = zeros(length(cSpeeds), 1);
for i=1:length(cSpeeds)
    slsParams.cSpeed_ = cSpeeds(i);
    
    if cSpeeds(i) == 1
        % in this case comm speed too slow; use approx instead of direct
        slsParams.mode_     = SLSMode.ApproxDLocalized;
        slsParams.robCoeff_ = 1000;
    else
        slsParams.mode_ = SLSMode.DLocalized;
    end
        
    slsOutsC        = state_fdbk_sls(sys, slsParams);
    clnorms(i)      = slsOutsC.clnorm_;
    
    if ismember(cSpeeds(i), cPrints)   
        [xC, uC] = simulate_system(sys, slsParams, slsOutsC, simParams);
        plot_heat_map(xC, sys.B2*uC, ['Comm Speed = ', num2str(cSpeeds(i))]);
    end
end

figure;
p4=plot(cSpeeds, clnorms, 'o-');
set(gca, 'xdir', 'reverse');
title([int2str(sys.Nx), ' Node Chain']);
xlabel('\alpha'); ylabel('Localized H_2-Norm Cost');

%% varying state spread (alpha)
alphas   = linspace(0, 0.8, 10);
printIdx = [1 5 10];

slsParams.mode_   = SLSMode.DLocalized;
slsParams.d_      = 6;
slsParams.cSpeed_ = 2;
clnorms           = zeros(length(alphas), 1);
specRadii         = zeros(length(alphas), 1);
for i=1:length(alphas)
    generate_dbl_stoch_chain(sys, rho, actDens, alphas(i)); % update A, B2
    % Nu changed, so have to change these as well
    sys.C1  = [speye(sys.Nx); sparse(sys.Nu, sys.Nx)];
    sys.D12 = [sparse(sys.Nx, sys.Nu); speye(sys.Nu)];
    
    slsOutsAlpha = state_fdbk_sls(sys, slsParams);
    clnorms(i)   = slsOutsAlpha.clnorm_;
    specRadii(i) = max(abs(eig(full(sys.A))));

    if ismember(i, printIdx)
        [xAlpha, uAlpha] = simulate_system(sys, slsParams, slsOutsAlpha, simParams);
        plot_heat_map(xAlpha, sys.B2*uAlpha, ['\alpha = ', num2str(alphas(i))])
    end
end

figure;
p5=plot(alphas, clnorms, 'o-');
title([int2str(sys. Nx), ' Node Chain']);
xlabel('\alpha'); ylabel('Localized H_2-Norm Cost');

figure;
p6=plot(alphas, specRadii, 'o-');
title([int2str(sys. Nx), ' Node Chain']);
xlabel('\alpha'); ylabel('Spectral Radius');