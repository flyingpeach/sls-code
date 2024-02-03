clear; close all; clc; 

% generate sys.A, sys.B2
Nx = 10; alpha = 0.2; rho = 1; actDens = 1;
sys = generate_dbl_stoch_chain(Nx, rho, actDens, alpha); 

sys.Nw  = sys.Nx; 
sys.Nz  = sys.Nu + sys.Nx;
sys.B1  = eye(sys.Nx); 
sys.C1  = [speye(sys.Nx); sparse(sys.Nu, sys.Nx)];
sys.D11 = sparse(sys.Nz, sys.Nw);
sys.D12 = [sparse(sys.Nx, sys.Nu); speye(sys.Nu)];
sys.sanity_check();

% simulation setup
simParams           = SimParams;
simParams.tSim_     = 25;
simParams.w_        = zeros(sys.Nx, simParams.tSim_); % disturbance

dist = [-1 1 1 -1 -1 -1 1 1 1 1 -1 -1 -1 -1 1 1 1 1 1 -1 -1 -1 -1 -1 -1];
simParams.w_(5, :) = dist;

%% (1) standard sls controller (for comparison)
slsParams    = SLSParams();
slsParams.T_ = 20;
slsParams.add_objective(SLSObjective.H2, 1);         % H2 objective
slsParams.add_constraint(SLSConstraint.Locality, 3); % locality = 3-hops

clMaps1  = state_fdbk_sls(sys, slsParams);
ctrller1 = Ctrller.ctrller_from_cl_maps(clMaps1);

[x1, u1] = simulate_state_fdbk(sys, ctrller1, simParams);
plot_heat_map(x1, sys.B2*u1, 'Standard SLS');

fprintf('\nLargest state magnitude: %.2f\n', max(abs(vec(x1))));

%% (2) saturation controller
satParams = SatParams();
satParams.noiseStd_ = 0.1;
satParams.eta_      = [0.9 1];
satParams.tau_      = -2; % no anti-windup is used since system is not stable
satParams.xMax_     = 1.01;
satParams.uMax_     = Inf;

% Solve for closed-loop maps for blending
[clMaps_b1, clMaps_b2] = two_blend_sls(sys, slsParams, satParams);
ctrller_b1 = Ctrller.ctrller_from_cl_maps(clMaps_b1);
ctrller_b2 = Ctrller.ctrller_from_cl_maps(clMaps_b2);

[x2, u2] = simulate_two_blend(sys, ctrller_b1, ctrller_b2, simParams, satParams);
plot_heat_map(x2, sys.B2*u2, 'Blended SLS');

fprintf('\nLargest state magnitude: %.2f\n', max(abs(vec(x2))));