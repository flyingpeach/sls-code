function plot_graph_animation(adjMtx, nodeCoords, x_, Bu_, waitTime, logScale)
% Plots topology of graph and animates states values at nodes
% Inputs
%   adjMtx     : adjacency matrix
%   nodeCoords : x,y coordinates of each node (in order)
%   x, Bu      : state and actuation values at nodes
%   waitTime   : amount of time to wait between steps
%   logScale   : whether to use the same scale as heat map plotter

if nargin == 5
    logScale = false;
end

minCoords = min(nodeCoords);
maxCoords = max(nodeCoords);
width     = maxCoords(1) - minCoords(1);
height    = maxCoords(2) - minCoords(2);
buffer    = 3;    % min size 3 cm
maxHeight = 15;   % max height 15 cm

if width > height
    height = height / width;
    width  = 1;
else
    width  = width / height;
    height = 1;
end

height = max(height * maxHeight, buffer); 
width  = 2 * max(width * maxHeight, buffer) + buffer; % two windows

% [left bottom width height]
figure('units','centimeters','outerposition',[1 1 width height])

Nx   = size(x_, 1);
TMax = size(x_, 2);

colormap default; 
cmap    = colormap; 

if logScale
    logmin  = -4; % these match the axes in plot_heat_map
    logmax  = 0;

    % clamp values within min and max
    normedx = max(min(log10(abs(x_)), logmax), logmin);
    normedu = max(min(log10(abs(Bu_)), logmax), logmin);
    
    % normalize values
    normedx = (normedx - logmin)./ (logmax - logmin);
    normedu = (normedu - logmin)./ (logmax - logmin);

else
    maxmagx = max(abs(vec(x_)));  % max magnitude
    maxmagu = max(abs(vec(Bu_)));
    normedx = abs(x_) ./ maxmagx;
    normedu = abs(Bu_) ./ maxmagu;
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