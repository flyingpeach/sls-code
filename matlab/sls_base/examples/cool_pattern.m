clear; close all; clc; 

% which things we want to plot
plotAnimation = false;
plotHeatMap   = true;

Nx = 50; alpha = 0.8; rho = 1; actDens = 1; 
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
simParams.tSim_     = 75;
simParams.w_        = zeros(sys.Nx, 100);
simParams.w_(floor(sys.Nx/2), 1) = 10;

% sls parameters
slsParams    = SLSParams();
slsParams.T_ = 9;

slsParams.add_objective(SLSObjective.H2, 1); 
slsParams.add_constraint(SLSConstraint.ActDelay, 1);
slsParams.add_constraint(SLSConstraint.CommSpeed, 1);
slsParams.add_constraint(SLSConstraint.Locality, 8);

slsParams.approx_      = true;
slsParams.stabCoeff_ = 1e3;

% find closed loop map + controller + simulate
clMaps  = state_fdbk_sls(sys, slsParams);
ctrller = Ctrller.ctrller_from_cl_maps(clMaps);
[x, u]  = simulate_state_fdbk(sys, ctrller, simParams);

nodeCoords = [zeros(1, sys.Nx);
              1:1:sys.Nx]';
          
onesVec = ones(sys.Nx - 1, 1);
adjMtx  = eye(sys.Nx) + diag(onesVec, 1) + diag(onesVec, -1);

if plotAnimation
    waitTime = 0.01;
    logScale = true;
    plot_graph_animation(adjMtx, nodeCoords, x, sys.B2*u, waitTime, logScale);
end

if plotHeatMap
    plot_heat_map(x, sys.B2*u, '');
end