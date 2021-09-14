function XiGCol = get_xig_col(G, Xi, s_cXi, cols, colG)
% cols are cols of Xi that the current processor is in charge of

gIdxs  = find(G(:, colG));

% sanity check that this processor has access to the desired cols in Xi
if ~ismember(gIdxs(1), cols) || ~ismember(gIdxs(2), cols)
    mpc_error('Tried to access a column of Xi we do not have access to');    
end

% Note: assuming uncoupled G, we are guaranteed that length(gIdxs) = 2
%       also we are guaranteed that the corresponding columns of Xi have
%       the same sparsity pattern, so we can directly stack them

XiCols = [Xi(s_cXi{gIdxs(1)}, gIdxs(1)) Xi(s_cXi{gIdxs(2)}, gIdxs(2))];
XiGCol = XiCols * G(gIdxs, colG);

end