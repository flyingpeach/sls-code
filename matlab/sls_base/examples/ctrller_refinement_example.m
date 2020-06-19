clear; close all; clc; 

% specify system matrices
sys    = LTISystem;
sys.Nx = 10; sys.Nw = sys.Nx; 

% generate sys.A, sys.B2
alpha = 0.4; rho = 1; actDens = 0.3;
generate_dbl_stoch_chain(sys, rho, actDens, alpha); 

sys.Nz  = sys.Nu + sys.Nx;
sys.B1  = eye(sys.Nx); 
sys.C1  = [speye(sys.Nx); sparse(sys.Nu, sys.Nx)];
sys.D11 = sparse(sys.Nz, sys.Nw);
sys.D12 = [sparse(sys.Nx, sys.Nu); speye(sys.Nu)];
sys.sanity_check();

slsParams       = SLSParams();
slsParams.T_    = 20;
slsParams.add_objective(SLSObjective.H2, 1);

% infinite horizon LQR
[P,L,K]   = dare(full(sys.A), full(sys.B2), eye(sys.Nx), eye(sys.Nu));
infLQRCost = 0;
for i=1:sys.Nx % sum of all costs from init condition
    x0    = zeros(sys.Nx, 1);
    x0(i) = 1;
    infLQRCost = infLQRCost + x0'*P*x0;
end

% centralized SLS
clMapsCent  = state_fdbk_sls(sys, slsParams);
ctrllerCent = Ctrller.ctrller_from_cl_maps(clMapsCent);

% distributed SLS (may be infeasible)
slsParams.add_constraint(SLSConstraint.CommSpeed, 1);
slsParams.add_constraint(SLSConstraint.Locality, 3);
clMapsDistr  = state_fdbk_sls(sys, slsParams);
ctrllerDistr = Ctrller.ctrller_from_cl_maps(clMapsDistr);

% approximate SLS (i.e. virtually localized)
slsParams.approx_    = true;
slsParams.stabCoeff_ = 1e3;
clMapsVirt  = state_fdbk_sls(sys, slsParams);
ctrllerVirt = Ctrller.ctrller_from_cl_maps(clMapsVirt);

% two-step method with controller L1 objective
slsParams.objectives_  = {}; % reset objectives

slsParams.stabCoeff_ = 1; % stability objective
eqnErrCoeff          = 1; % equation error penalty
slsParams.add_objective(SLSObjective.L1, 1); % controller L1 objective
ctrller2Step_1 = refine_ctrller(sys, clMapsCent, slsParams, eqnErrCoeff); 

% lower order two-step method
slsParams.T_ = 2;
ctrller2Step_2 = refine_ctrller(sys, clMapsCent, slsParams, eqnErrCoeff);

%% calculate + print stats for comparison
tTotal = 200;

[lqrCostCent, intRadCent, l1normCent] = get_ctrller_stats(sys, ctrllerCent, tTotal);
fprintf('Centralized    : cost=%.3f, rad=%.3f, l1norm=%.3f\n', ...
        lqrCostCent / infLQRCost, intRadCent, l1normCent);

if strcmp(clMapsDistr.solveStatus_, 'Solved')
    [lqrCostDistr, intRadDistr, l1normDistr] = get_ctrller_stats(sys, ctrllerDistr, tTotal); 
    fprintf('Distributed    : cost=%.3f, rad=%.3f, l1norm=%.3f\n', ...
            lqrCostDistr / infLQRCost, intRadDistr, l1normDistr);
else
    fprintf('Distributed    : Infeasible!\n');
end

[lqrCostVirt, intRadVirt, l1normVirt] = get_ctrller_stats(sys, ctrllerVirt, tTotal);
fprintf('Virtualized    : cost=%.3f, rad=%.3f, l1norm=%.3f\n', ...
        lqrCostVirt / infLQRCost, intRadVirt, l1normVirt);

[lqrCost2Step_1, intRad2Step_1, l1norm2Step_1] = get_ctrller_stats(sys, ctrller2Step_1, tTotal);
fprintf('Two-step (T=20): cost=%.3f, rad=%.3f, l1norm=%.3f\n', ...
         lqrCost2Step_1 / infLQRCost, intRad2Step_1, l1norm2Step_1);

[lqrCost2Step_2, intRad2Step_2, l1norm2Step_2] = get_ctrller_stats(sys, ctrller2Step_2, tTotal);
fprintf('Two-step (T=2) : cost=%.3f, rad=%.3f, l1norm=%.3f\n', ...
         lqrCost2Step_2 / infLQRCost, intRad2Step_2, l1norm2Step_2);
