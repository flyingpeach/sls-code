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

plot_graph(adjMtx, nodeCoords, 'w');

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
comms   = 3;  % communication speed
ta      = 1;  % actuation delay
TFIR    = 10; % finite impulse response horizon
TMax    = 20; % amount of time to simulate
TSim    = TFIR + TMax;

[R, M]  = sf_sls_d_localized(A, B, C, D, TFIR, d, comms, ta, 'H2');

loc = 4; % where disturbance hits

% useful for visualizing time; not for visualizing locality
% make_heat_map(A,B,TFIR,Nx,Nu,R,M,loc,TSim,'Localized')

% TODO: code overlap
% this code is copied from make_heat_map with some modifications
% it simulates the system for TSim time

w_d    = zeros(Nx,TSim);
Tstart = TFIR+1;
B1     = eye(Nx);

w_d(loc,Tstart) = 10;

x     = zeros(Nx,TSim);
u     = zeros(Nu,TSim);
x_ref = zeros(Nx,TSim);
w_est = zeros(Nx,TSim);

MATx = []; MATu = [];

for i=Tstart:1:TSim-1

    w_est(:,i-1) = x(:,i) - x_ref(:,i);
    
    for jj=1:1:TFIR
        u(:,i) = u(:,i) + M{jj}*w_est(:,i-jj);
    end
    
    x(:,i+1) = A*x(:,i) + B1*w_d(:,i)+ B*u(:,i);
    
    for jj=2:1:TFIR
        x_ref(:,i+1) = x_ref(:,i+1) + R{jj}*w_est(:,i+1-jj);
    end

    MATx = [MATx,log10(abs(x(:,i)))];
    MATu = [MATu,log10(abs(B*u(:,i)))];   
end

% visualize disturbance propagation
colormap jet; 
cmap    = colormap; 

finMATx = MATx(isinf(MATx) == 0); % vector with non-infinite values of x
finMATu = MATu(isinf(MATu) == 0); 

minx = min(finMATx); % min value
minu = min(finMATu);

maxmagx = max(abs(finMATx)); % max magnitude
maxmagu = max(abs(finMATu)); % max magnitude

subplot(2,1,1);
plot_graph(adjMtx, nodeCoords, cmap(1,:));
title('normalized log(x)')
colorbar;

subplot(2,1,2);
plot_graph(adjMtx, nodeCoords, cmap(1,:));
colorbar;
title('normalized log(u)')

for time=1:TMax-1
    pause(0.5);
    if time > 1
        delete(timeText)
    end
    normedx = (MATx(:,time) - minx) ./ maxmagx;
    normedu = (MATu(:,time) - minu) ./ maxmagu;

    if (time == TFIR)
        timeText = text(2, -0.3, strcat('t=', num2str(time)), 'Color', 'r');
    else
        timeText = text(2, -0.3, strcat('t=', num2str(time)));
    end
        
    for node=1:Nx
        subplot(2,1,1)
        plot_vertex(node, nodeCoords, get_colour(normedx(node), cmap));
        subplot(2,1,2)
        plot_vertex(node, nodeCoords, get_colour(normedu(node), cmap));
    end
end
