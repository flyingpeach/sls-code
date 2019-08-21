function [colour] = get_colour(val, cmap)
% maps value to a colour in the colour map
%   val    : normalized value (or inf)
%   cmap   : colormap
%   colour : colour as RGB coordinate

cmapval = linspace(0, 1, length(cmap));

if isinf(val)
    if sign(val) == -1
        idx = 1; % set to lowest
    else
        idx = length(cmap); % set to highest
    end
else
    [junk, idx] = min(abs(cmapval - val));
end

colour = cmap(idx,:);