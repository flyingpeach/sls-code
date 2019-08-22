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

figure(1)
plot_graph(adjMtx, nodeCoords, 'w');

%% control / simulation
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

    MATx = [MATx, x(:,i)];
    MATu = [MATu, B*u(:,i)];   
end

%% visualization 1: heat map
% useful for visualizing time; not for visualizing locality
figure(2)
make_heat_map(A,B,TFIR,Nx,Nu,R,M,loc,TSim,'Localized')

%% visualization 2: graph-based animation

figure(3)

colormap jet; 
cmap    = colormap; 

maxmagx = max(abs(vec(MATx))); % max magnitude
maxmagu = max(abs(vec(MATu))); % max magnitude

subplot(2,1,1);
plot_graph(adjMtx, nodeCoords, cmap(1,:));
title('normalized x')
colorbar;

subplot(2,1,2);
plot_graph(adjMtx, nodeCoords, cmap(1,:));
colorbar;
title('normalized u')

for time=1:TMax-1
    pause(0.5);
    if time > 1
        delete(timeText)
    end
    normedx = abs(MATx(:,time)) ./ maxmagx;
    normedu = abs(MATu(:,time)) ./ maxmagu;

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

%% visualization 3: time-plots per node

maxy = max(max(vec(MATx)), max(vec(MATu))) + 2;
miny = min(min(vec(MATx)), min(vec(MATu))) - 2;

figure(4)
for node=1:Nx
    subplot(Nx, 1, node);
    hold on
    stairs([1:TMax-1], MATx(node,:));
    stairs([1:TMax-1], MATu(node,:));
    set(gca,'XTickLabel',[])
    ylabel(num2str(node));
    ylim([miny maxy]);
end

legend('x','u');
xlabel('time step');
set(gca,'XTickLabelMode','auto');

