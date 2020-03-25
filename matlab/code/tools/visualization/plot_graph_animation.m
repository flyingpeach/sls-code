function plot_graph_animation(adjMtx, nodeCoords, slsParams, x, Bu, waitTime, logScale)
% Plots topology of graph and animates states values at nodes
% Inputs
%   adjMtx     : adjacency matrix
%   nodeCoords : x,y coordinates of each node (in order)
%   slsParams  : SLSParams containing parameters
%   x, Bu      : state and actuation values at nodes
%   waitTime   : amount of time to wait between steps
%   logScale   : whether to use the same scale as heat map plotter

if nargin == 6
    logScale = false;
end

% make full-screen
figure('units','normalized','outerposition',[0 0.1 1 0.7])

Nx   = size(x, 1);
TMax = size(x, 2);

colormap default; 
cmap    = colormap; 

if logScale
    logmin  = -4; % these match the axes in plot_heat_map
    logmax  = 0;

    % clamp values within min and max
    normedx = max(min(log10(abs(x)), logmax), logmin);
    normedu = max(min(log10(abs(Bu)), logmax), logmin);
    
    % normalize values
    normedx = (normedx - logmin)./ (logmax - logmin);
    normedu = (normedu - logmin)./ (logmax - logmin);

else
    maxmagx = max(abs(vec(x)));  % max magnitude
    maxmagu = max(abs(vec(Bu)));
    normedx = abs(x) ./ maxmagx;
    normedu = abs(Bu) ./ maxmagu;
end

subplot(1,2,1);
plot_graph(adjMtx, nodeCoords, cmap(1,:));
title('state (normalized)')
colorbar;

subplot(1,2,2);
plot_graph(adjMtx, nodeCoords, cmap(1,:));
colorbar;
title('control (normalized)')

for time=1:TMax-1
    pause(waitTime);

    for node=1:Nx
        subplot(1,2,1)
        plot_vertex(node, nodeCoords, get_colour(normedx(node, time), cmap));
        subplot(1,2,2)
        plot_vertex(node, nodeCoords, get_colour(normedu(node, time), cmap));
    end
end
end


% local functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [colour] = get_colour(val, cmap)
% Maps value to a colour in the colour map
%   val    : normalized value (or inf)
%   cmap   : colormap
%   colour : colour as RGB coordinate

cmapDim = size(cmap);

cmapval = linspace(0, 1, cmapDim(1));

if isinf(val)
    if sign(val) == -1
        idx = 1; % set to lowest
    else
        idx = cmapDim(1); % set to highest
    end
else
    [junk, idx] = min(abs(cmapval - val));
end

colour = cmap(idx,:);
end