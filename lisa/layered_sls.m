clear;
clc;

% which things we want to plot
plot_animation = false;
plot_time_traj = true;

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
    fprintf('[error] adjacency matrix is not symmetric!\n')
end

% dynamics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx = length(adjMtx); % number of nodes (and therefore states)
Nu = Nx; % for now fully actuated

% xi[t+1] = aii*xi[t] + sum(aij*xj[t]) + bii*ui[t] + di
% where xj are neighbours of xi
aii = 1;      % self effect
bii = 1;      % control effect
aij = 2 / Nx; % neighbour effect

A = aii * eye(Nx) + aij * adjMtx;
B = bii * eye(Nx);

% sls setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
obj = 'traj_track';

C = [speye(Nx); sparse(Nu,Nx)];
D = [sparse(Nx,Nu); speye(Nu)];

d       = 3;  % d-hop locality constraint
comms   = 3;  % communication speed
ta      = 1;  % actuation delay
TFIR    = 17; % finite impulse response horizon
TMax    = 25; % amount of time to simulate

% desired trajectory %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xDes = zeros(Nx, TFIR-2);

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

% pad first column with zeros
% also pad with zeros if xDes not specified for all time
sizexd   = size(xDes);
timexd   = sizexd(2);
timediff = TMax - timexd - 1;
xDes     = [zeros(Nx, 1) xDes zeros(Nx, timediff-1)];

% sls and simulate system %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[R, M] = sf_sls_basic(A, B, C, D, TFIR, obj, xDes);
%[R, M]  = sf_sls_d_localized(A, B, C, D, TFIR, d, comms, ta, obj, xDesired);

% TODO: code overlap with make_heat_map
% simulate system as per equation (2.8)

w       = zeros(Nx, TMax); % true disturbance
w(6, 1) = 10;

B1  = eye(Nx); % see eqn (3.1);  transfer fn from noise to state

x     = zeros(Nx,TMax); u     = zeros(Nu,TMax);
x_hat = zeros(Nx,TMax); w_hat = zeros(Nx,TMax);

for t=1:1:TMax-1
    for tau=1:1:min(t-1, TFIR)
       u(:,t) = u(:,t) + M{tau}*w_hat(:,t-tau);
    end

    for tau=1:1:min(t-1, TFIR-1)
       x_hat(:,t+1) = x_hat(:,t+1) + R{tau+1}*w_hat(:,t-tau);       
    end    

    x(:,t+1) = A*x(:,t) + B1*w(:,t)+ B*u(:,t);
    w_hat(:,t) = x(:,t+1) - x_hat(:,t+1);
end

u = B*u; % want to look at actuation at each node, not the u themselves

% visualization: graph-based animation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_animation
    figure(1)
    colormap jet; 
    cmap    = colormap; 

    maxmagx = max(abs(vec(x))); % max magnitude
    maxmagu = max(abs(vec(u))); % max magnitude

    subplot(2,1,1);
    plot_graph(adjMtx, nodeCoords, cmap(1,:));
    title('normalized x')
    colorbar;

    subplot(2,1,2);
    plot_graph(adjMtx, nodeCoords, cmap(1,:));
    colorbar;
    title('normalized u')

    for time=1:TMax-1
        pause(0.2);
        if time > 1
            delete(timeText)
        end
        normedx = abs(x(:,time)) ./ maxmagx;
        normedu = abs(u(:,time)) ./ maxmagu;

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
end

% visualization: time trajectory plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if plot_time_traj
    figure(2)
    % TODO: hack: need to shift xDesired
    xDes = [zeros(Nx, 1) xDes(:, 1:timexd + timediff)];
    
    maxy = max([max(vec(x)) max(vec(u)) max(vec(xDes))]) + 2;
    miny = min([min(vec(x)) min(vec(u)) max(vec(xDes))]) - 2;

    err = abs(xDes - x) / max(vec(xDes));
    maxe = max(vec(err)) * 1.1; 
    mine = min(vec(err));
    
    for node=1:Nx
        subplot(Nx, 2, node * 2 - 1);
        hold on
        stairs(1:TMax, x(node,:));
        stairs(1:TMax, xDes(node,:));
        stairs(1:TMax, u(node,:));
        set(gca,'XTickLabel',[]);
        ylabel(num2str(node));
        ylim([miny maxy]);
        
        subplot(Nx, 2, node * 2);
        stairs(1:TMax, err(node,:));
        set(gca,'XTickLabel',[]);
        ylabel(num2str(node));
        ylim([mine maxe]);
    end

    subplot(Nx, 2, Nx * 2 - 1)
    legend('x', 'xDes', 'u');
    xlabel('time step');
    set(gca,'XTickLabelMode','auto');

    subplot(Nx, 2, Nx * 2)
    legend('normed err');
    xlabel('time step');
    set(gca,'XTickLabelMode','auto');

 end