function x = col_robust_closed(c, m, g, v)
% The arguments are consistent with the proximal operator solution
% This general expression applies to both column-partitions
% of the column-wise update

x = c*[m'*v; g];
end