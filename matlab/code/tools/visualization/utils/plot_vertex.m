function plot_vertex(node, nodeCoords, colour)
% Adds colour & label to specified vertex / node
%   node       : the node we're plotting
%   nodeCoords : x,y coordinates of each node (in order)
%   colour     : colour of node, either a letter or a RGB coordinate

% make node a different colour
plot(nodeCoords(node,1), nodeCoords(node,2), ...
     'o','MarkerSize', 10, 'MarkerFaceColor',colour,'MarkerEdgeColor','k');