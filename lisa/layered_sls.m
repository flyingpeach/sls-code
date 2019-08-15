%% parameters
Nx = 7; % number of nodes (and therefore states)
Nu = 7; % for now fully actuated

% xi[t+1] = aii*xi[t] + sum(aij*xj[t]) + bii*ui[t] + di
% where xj are neighbours of xi
aii = 1;      % self effect
bii = 1;      % control effect
aij = 2 / Nx; % neighbour effect

%% specify and visualize graph structure
% x,y coordinates of each node (in order)
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

% dynamic axis limits
minCoords = min(nodeCoords);     maxCoords = max(nodeCoords); 
xlowerlim = minCoords(1) - 0.5;  xupperlim = maxCoords(1) + 0.5;
ylowerlim = minCoords(2) - 0.5;  yupperlim = maxCoords(2) + 0.5;

gplot(adjMtx, nodeCoords, '-ok');
hold on
for i=1:Nx % make nodes a different colour
    plot(nodeCoords(i,1), nodeCoords(i,2), ...
         'o','MarkerFaceColor','g','MarkerEdgeColor','k');
end
axis([xlowerlim xupperlim ylowerlim yupperlim])
hold off

%% dynamics
A = aii * eye(Nx) + aij * adjMtx;
B = bii * eye(Nx);