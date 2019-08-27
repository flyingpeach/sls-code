function plot_graph_animation(adjMtx, nodeCoords, params, x, Bu)
% Plots topology of graph and animates states values at nodes
% Inputs
%   adjMtx     : adjacency matrix
%   nodeCoords : x,y coordinates of each node (in order)
%   params     : SLSParams containing parameters
%   x, Bu      : state and actuation values at nodes

figure;
Nx   = size(x, 1);
TMax = size(x, 2);

colormap jet; 
cmap    = colormap; 

maxmagx = max(abs(vec(x)));  % max magnitude
maxmagu = max(abs(vec(Bu)));

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
    normedu = abs(Bu(:,time)) ./ maxmagu;

    if (time == params.tFIR_)
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

% local function
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