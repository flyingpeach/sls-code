function plot_heat_map(x, Bu, myTitle)
% Plots log-based heat map for x, u
% Inputs
%    x, Bu    : state and actuation values at nodes
%    myTitle  : overall title of the heat maps

figure; suptitle(myTitle);

if Bu == zeros(size(Bu, 1), size(Bu, 2));
    % don't subplot; plot only x
    imagesc(log10(abs(x))); colorbar; colormap jet; caxis([-4 0]);
    title('log10(|x|)'); xlabel('Time'); ylabel('Space');
else
    subplot(1,2,1)
    imagesc(log10(abs(x))); colorbar; colormap jet; caxis([-4 0]);
    title('log10(|x|)'); xlabel('Time'); ylabel('Space');

    subplot(1,2,2)
    imagesc(log10(abs(Bu))); colorbar; colormap jet; caxis([-4 0]);
    title('log10(|u|)'); xlabel('Time');
end

% These are the font settings from previous papers; uncomment as wanted
% set(gca,'FontSize',16,'fontWeight','bold')
% set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','bold')

end

