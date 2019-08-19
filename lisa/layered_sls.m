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

%% control parameters
Nx = length(adjMtx); % number of nodes (and therefore states)
Nu = Nx; % for now fully actuated

% xi[t+1] = aii*xi[t] + sum(aij*xj[t]) + bii*ui[t] + di
% where xj are neighbours of xi
aii = 1;      % self effect
bii = 1;      % control effect
aij = 2 / Nx; % neighbour effect


A = aii * eye(Nx) + aij * adjMtx;
B = bii * eye(Nx);


