function heat_map_plotter(mtx, myTitle, myXLabel, myYLabel)
% plots the individual heatmaps using the 
%   mtx     : the matrix containing the info to be plotted

imagesc((mtx));
colorbar; colormap jet; caxis([-4 0]);
title(myTitle); xlabel(myXLabel); ylabel(myYLabel);
set(gca,'FontSize',16,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','bold')

end