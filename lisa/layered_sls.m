%% specify and visualize graph structure
clear;
clc;
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
          
% sanity check
if ~issymmetric(adjMtx)
    fprintf('[error] adjacency matrix is not symmetric!\n')
end

plot_graph(adjMtx, nodeCoords);

%% control
Nx = length(adjMtx); % number of nodes (and therefore states)
Nu = Nx; % for now fully actuated

% xi[t+1] = aii*xi[t] + sum(aij*xj[t]) + bii*ui[t] + di
% where xj are neighbours of xi
aii = 1;      % self effect
bii = 1;      % control effect
aij = 2 / Nx; % neighbour effect

A = aii * eye(Nx) + aij * adjMtx;
B = bii * eye(Nx);

C = [speye(Nx); sparse(Nu,Nx)];
D = [sparse(Nx,Nu); speye(Nu)];

d       = 3;  % d-hop locality constraint
comms   = 20; % computation speed
ta      = 0;  % actuation delay
TFIR    = 20; % finite impulse response horizon
TAfter  = 20; % amount of time to simulate after FIR horizon
TSim    = TFIR + TAfter;

[R, M]  = sf_sls_d_localized(A, B, C, D, TFIR, d, comms, ta, 'H2');

loc = floor(Nx/2); % where disturbance hits

% note: right now this is only useful to visualize FIR (time)
% locality is not reflected
make_heat_map(A, B, TFIR, Nx, Nu, R, M, loc, TSim, 'Localized')

