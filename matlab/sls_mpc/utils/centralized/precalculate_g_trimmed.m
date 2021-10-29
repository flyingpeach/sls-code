function GTrims = precalculate_g_trimmed(G, s_rXi, nH)
% Fetch zero-trimmed submatrices of G (to be postmultipled by Xi) 
% Indexed by rows of Xi/H

GTrims = cell(nH, 1);
for rowH=1:nH
   G_          = G(s_rXi{rowH}, :);
   zeroCols    = find(all(G_ == 0, 1));
   keepCols    = setdiff(1:size(G_,2), zeroCols);      
   GTrims{rowH} = G_(:, keepCols);
end