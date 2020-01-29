clear; close all; clc; 

% which things we want to plot
plotAnimation = true;
plotTimeTraj  = false;
plotHeatMap   = false;

% graph architecture %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nodeCoords = [0 1;
              1 1;
              3 1;
              4 1;
              0 0;
              3 0;
              4 0;];

adjMtx = [0 1 0 0 1 0 0;
          1 0 0 0 1 0 1;
          0 0 0 1 1 1 1;
          0 0 1 0 0 1 1;
          1 1 1 0 0 0 0;
          0 0 1 1 0 0 1;
          0 1 1 1 0 1 0];
          
% sanity check on graph architecture
if ~issymmetric(adjMtx)
    error('Adjacency matrix is not symmetric!')
end

% dynamics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xi[t+1] = aii*xi[t] + sum(aij*xj[t]) + bii*ui[t] + di
% where xj are neighbours of xi
aii = 1;      % self effect
bii = 1;      % control effect
aij = 2 / size(adjMtx, 1); % neighbour effect

% create a system
sys     = LTISystem;
sys.Nx  = size(adjMtx, 1); % one state per node
sys.Nu  = sys.Nx;          % fully actuated
sys.A   = aii * eye(sys.Nx) + aij * adjMtx;
sys.B1  = eye(sys.Nx);
sys.B2  = bii * eye(sys.Nx);
sys.C1  = [speye(sys.Nx); sparse(sys.Nu, sys.Nx)];
sys.D12 = [sparse(sys.Nx, sys.Nu); speye(sys.Nu)];

% sls %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
slsParams           = SLSParams;
slsParams.actDelay_ = 1;
slsParams.cSpeed_   = 3;
slsParams.d_        = 3;
slsParams.T_     = 17;
slsParams.obj_      = Objective.H2;
slsParams.mode_     = SLSMode.Basic;
slsParams.rfd_      = false;

% note: we are doing H2-SLS on the error of the trajectory (xd, ud)
slsOuts = state_fdbk_sls(sys, slsParams);

% desired trajectory %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xDes    = zeros(sys.Nx, slsParams.T_);

xDes(6:7, 1)  = 1; % pedal
xDes(6:7, 5)  = 1;
xDes(6:7, 9)  = 1;
xDes(6:7, 13) = 1;

xDes(3, 1)  = 1; % melody
xDes(3, 2)  = 1;
xDes(4, 3)  = 1;
xDes(5, 4)  = 1;
xDes(5, 5)  = 1;
xDes(4, 6)  = 1;
xDes(3, 7)  = 1;
xDes(2, 8)  = 1;
xDes(1, 9)  = 1;
xDes(1, 10) = 1;
xDes(2, 11) = 1;
xDes(3, 12) = 1;
xDes(3, 13) = 1;
xDes(2, 14) = 1;
xDes(2, 15) = 1;

trajMag = 10; % magnitude of trajectory (const magnitude for now)
xDes    = xDes .* trajMag;

% use LQR to find xd, ud
xPenalty = 1;
uPenalty = 0.1;
QLQR     = xPenalty * eye(sys.Nx);
RLQR     = uPenalty * eye(sys.Nu);

cvx_begin quiet
variable xd(sys.Nx, slsParams.T_)
variable ud(sys.Nu, slsParams.T_)

objective = 0;
for t = 1:slsParams.T_
    objective = objective + norm(QLQR*(xd(:,t)-xDes(:,t))) + norm(RLQR*ud(:,t));
    if t < slsParams.T_
        xd(:,t+1) == sys.A*xd(:,t) + sys.B2*ud(:,t); % enforce dynamics
    end
end

minimize(objective);
cvx_end

lqrStatus = sprintf('LQR trajectory tracking, xPenalty=%0.2f, uPenalty=%0.2f cost=%0.2f', ...
                    xPenalty, uPenalty, objective);
disp(lqrStatus);

% simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simParams           = SimParams;
simParams.tSim_     = 25;
simParams.openLoop_ = false;
% white noise disturbance
simParams.w_        = wgn(sys.Nx, simParams.tSim_, trajMag ./ 100);

% note: these are the errors of the trajectories
[xerr, uerr]  = simulate_system(sys, simParams, slsOuts.R_, slsOuts.M_);

% pad trajectories with zeros so that dimensions work out
timediff = simParams.tSim_ - slsParams.T_;
x        = xerr + [xd zeros(sys.Nx, timediff)]; % x = xerror + trajectory
u        = uerr + [ud zeros(sys.Nu, timediff)]; % u = uerror + trajectory
xDes     = [xDes zeros(sys.Nx, timediff)];

% visualizations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Bu = sys.B2*u; % want to look at actuation at each node

if plotAnimation
    waitTime = 0.5;
    plot_graph_animation(adjMtx, nodeCoords, slsParams, x, Bu, waitTime);
end

if plotTimeTraj
   plot_time_traj(x, Bu, xDes);
end
 
if plotHeatMap
    plot_heat_map(x, Bu, '');
end
 