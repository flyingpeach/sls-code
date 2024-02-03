function plot_heat_map(x, Bu, myTitle)
% Plots log-based heat map for x, u
% Inputs
%    x, Bu    : state and actuation values at nodes
%    myTitle  : overall title of the heat maps

figure;
logmin = -4;
logmax = 0;

subplot(1,2,1)
imagesc(log10(abs(x))); colorbar; colormap jet; caxis([logmin logmax]);
title([myTitle ' log10(|x|)']); xlabel('Time'); ylabel('Space');
pos = get(gca, 'Position');
set(gca, 'Position', [pos(1) pos(2)+0.01 pos(3) pos(4)-0.01]);

subplot(1,2,2)
imagesc(log10(abs(Bu))); colorbar; colormap jet; caxis([logmin logmax]);
title([myTitle ' log10(|u|)']); xlabel('Time');

% These are the font settings from previous papers; uncomment as wanted
% set(gca,'FontSize',16,'fontWeight','bold')
% set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','bold')

end

