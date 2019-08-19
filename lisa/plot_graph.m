function plot_graph(adjMtx, nodeCoords)
% visualizes graph
%   adjMtx    : adjacency matrix
%   nodeCoords: x,y coordinates of each node (in order)

% dynamic axis limits
minCoords = min(nodeCoords);     maxCoords = max(nodeCoords); 
xlowerlim = minCoords(1) - 0.5;  xupperlim = maxCoords(1) + 0.5;
ylowerlim = minCoords(2) - 0.5;  yupperlim = maxCoords(2) + 0.5;

gplot(adjMtx, nodeCoords, '-ok');

hold on
for i=1:length(adjMtx) % make nodes a different colour
    plot(nodeCoords(i,1), nodeCoords(i,2), ...
         'o','MarkerFaceColor','g','MarkerEdgeColor','k');
end
axis([xlowerlim xupperlim ylowerlim yupperlim])
hold off