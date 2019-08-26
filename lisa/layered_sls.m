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
slsparams           = SLSParams;
slsparams.actDelay_ = 1;
slsparams.cSpeed_   = 3;
slsparams.d_        = 3;
slsparams.tFIR_     = 17;
slsparams.obj_      = Objective.TrajTrack;

slsparams.setDesiredTraj(xDes);

% sls and simulate system %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[R, M] = sf_sls_basic(sys, slsparams);

% amount of time to simulate
TMax    = 25; 
% disturbance
w       = zeros(sys.Nx, TMax);
w(6, 1) = 10;

[x, u]  = simulate_system(sys, slsparams, TMax, R, M, w);

u = sys.B2*u; % want to look at actuation at each node, not the u themselves

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

        for node=1:sys.Nx
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
    % TODO: hack: need to shift xDesired twice since we already shifted it
    %             once in SLSParams
    timediff = TMax - size(xDes, 2);
    xDes     = [zeros(sys.Nx, 2) xDes zeros(sys.Nx, timediff-2)];
    
    maxy = max([max(vec(x)) max(vec(u)) max(vec(xDes))]) + 2;
    miny = min([min(vec(x)) min(vec(u)) max(vec(xDes))]) - 2;

    err = abs(xDes - x) / max(vec(xDes));
    maxe = max(vec(err)) * 1.1; 
    mine = min(vec(err));
    
    for node=1:sys.Nx
        subplot(sys.Nx, 2, node * 2 - 1);
        hold on
        stairs(1:TMax, x(node,:));
        stairs(1:TMax, xDes(node,:));
        stairs(1:TMax, u(node,:));
        set(gca,'XTickLabel',[]);
        ylabel(num2str(node));
        ylim([miny maxy]);
        
        subplot(sys.Nx, 2, node * 2);
        stairs(1:TMax, err(node,:));
        set(gca,'XTickLabel',[]);
        ylabel(num2str(node));
        ylim([mine maxe]);
    end

    subplot(sys.Nx, 2, sys.Nx * 2 - 1)
    legend('x', 'xDes', 'u');
    xlabel('time step');
    set(gca,'XTickLabelMode','auto');

    subplot(sys.Nx, 2, sys.Nx * 2)
    legend('normed err');
    xlabel('time step');
    set(gca,'XTickLabelMode','auto');

 end