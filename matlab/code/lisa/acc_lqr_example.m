clear; close all; clc; 

% set up system
sys    = LTISystem();
sys.Nx = 10;

alpha = 0.4; rho = 1; actDens = 0.3;
generate_dbl_stoch_chain(sys, rho, actDens, alpha); % generate sys.A, sys.B2

sys.B1  = eye(sys.Nx);
sys.C1  = [speye(sys.Nx); sparse(sys.Nu, sys.Nx)];
sys.D12 = [sparse(sys.Nx, sys.Nu); speye(sys.Nu)];

% sls parameters
slsParams           = SLSParams();
slsParams.tFIR_     = 20;
slsParams.obj_      = Objective.H2;
slsParams.mode_     = SLSMode.Basic;

slsOutsCent = state_fdbk_sls(sys, slsParams);

% infinite horizon
[P,L,K] = dare(full(sys.A), full(sys.B2), eye(sys.Nx), eye(sys.Nu));

infH2Cost = 0;
for i=1:sys.Nx % H2 cost is sum of all costs from init condition
    x0    = zeros(sys.Nx, 1);
    x0(i) = 1;
    infH2Cost = infH2Cost + x0'*P*x0;
end

% approx d-localized SLS
slsParams.mode_     = SLSMode.Localized;
slsParams.d_        = 3;
slsParams.robCoeff_ = 1e3;
slsParams.approx_   = true;

tTotal = slsParams.tFIR_+100;

slsOutsApproxLoc = state_fdbk_sls(sys, slsParams);
ctrllerALoc      = Ctrller();
ctrllerALoc.Rc_  = slsOutsApproxLoc.R_;
ctrllerALoc.Mc_  = slsOutsApproxLoc.M_;

% our method for d-localized SLS
eps_base   = 2.22e-16;
eps_nullsp = eps_base.^(3/8);

cParams             = CtrllerParams();
cParams.mode_       = SLSMode.Localized;
cParams.eps_nullsp_ = eps_nullsp;
cParams.tFIR_       = slsParams.tFIR_;
cParams.d_          = slsParams.d_;
cParams.CLDiffPen_  = 1e2;

ctrllerOurs = find_ctrller(sys, slsParams, slsOutsCent, cParams);

cStatsALoc = CtrllerStats(eps_nullsp, cParams.tFIR_);
csALoc{1}  = ctrllerALoc;
cStatsALoc = get_ctrller_stats(cStatsALoc, sys, slsParams, slsOutsCent, csALoc);

cStatsOurs = CtrllerStats(eps_nullsp, cParams.tFIR_);
csOurs{1}  = ctrllerOurs;
cStatsOurs = get_ctrller_stats(cStatsOurs, sys, slsParams, slsOutsCent, csOurs);

disp(sprintf('Centralized : cost=%.3f, rad=%.3f, l1norm=%.3f', ...
             cStatsOurs.LQRCostOrig / infH2Cost, ...
             cStatsOurs.IntSpecRadiusOrig, ...
             cStatsOurs.L1NormOrig))

disp(sprintf('Approx Loc. : cost=%.3f, rad=%.3f, l1norm=%.3f', ...
             cStatsALoc.LQRCosts / infH2Cost, ...
             cStatsALoc.IntSpecRadii_c, ...
             cStatsALoc.L1Norms))

disp(sprintf('Our method  : cost=%.3f, rad=%.3f, l1norm=%.3f', ...
             cStatsOurs.LQRCosts / infH2Cost, ...
             cStatsOurs.IntSpecRadii_c, ...
             cStatsOurs.L1Norms))
