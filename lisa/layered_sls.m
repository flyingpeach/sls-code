clear; close all; clc; 

% which things we want to plot
plotAnimation = false;
plotTimeTraj  = true;
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

% desired trajectory %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xDes = zeros(sys.Nx, 15);

xDes(6:7, 1)  = 20; % pedal
xDes(6:7, 5)  = 20;
xDes(6:7, 9)  = 20;
xDes(6:7, 13) = 20;

xDes(3, 1)  = 20; % melody
xDes(3, 2)  = 20;
xDes(4, 3)  = 20;
xDes(5, 4)  = 20;
xDes(5, 5)  = 20;
xDes(4, 6)  = 20;
xDes(3, 7)  = 20;
xDes(2, 8)  = 20;
xDes(1, 9)  = 20;
xDes(1, 10) = 20;
xDes(2, 11) = 20;
xDes(3, 12) = 20;
xDes(3, 13) = 20;
xDes(2, 14) = 20;
xDes(2, 15) = 20;

% sls setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
slsParams           = SLSParams;
slsParams.actDelay_ = 1;
slsParams.cSpeed_   = 3;
slsParams.d_        = 3;
slsParams.tFIR_     = 17;
slsParams.obj_      = Objective.TrajTrack;
slsParams.mode_     = SLSMode.Basic;
slsParams.rfd_      = false;

slsParams.setDesiredTraj(xDes);

% sls and simulate system %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation params
TMax    = 25;                  % amount of time to simulate
w       = zeros(sys.Nx, TMax); % disturbance
w(3, 1) = 10;

slsOuts = state_fdbk_sls(sys, slsParams);
[x, u]  = simulate_system(sys, slsParams, slsOuts, TMax, w);

% visualizations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Bu = sys.B2*u; % want to look at actuation at each node

if plotAnimation
    plot_graph_animation(adjMtx, nodeCoords, slsParams, x, Bu);
end

if plotTimeTraj
   % TODO: hack: need to shift xDesired twice since we already shifted it
   %             once in SLSParams
   timediff = TMax - size(xDes, 2);
   xDes     = [zeros(sys.Nx, 2) xDes zeros(sys.Nx, timediff-2)];
   plot_time_traj(x, Bu, xDes);
end
 
if plotHeatMap
    plot_heat_map(x, Bu, '');
end
 