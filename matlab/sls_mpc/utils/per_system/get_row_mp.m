function [MPRow, s_r] = get_row_mp(M, P, s_rP, rows, rows_neighbor, row)
% Returns row of matrix M*Psi or M*Phi (generally M is H or C)
% rows: rows of phi/psi that the current processor is in charge of
% P   : can be either Phi or Psi
% s_rP: can be s_rPhi or s_rPsi (corresponding to P)
% rows_neighbor: rows of psi that other processors are in charge of
%                that the current processor can access (i.e. neighbors)

pRows = find(M(row, :)); % all rows accessed via H*P
 
% sanity check that this processor has access to the desired row in Psi
if ~all(ismember(pRows, [rows, rows_neighbor]))
    mpc_error('Tried to access a row of Phi/Psi we do not have access to');
end

% Different rows will have different s_rPsi, take the union
s_r = [];

for pRow = pRows
    s_r = union(s_r, s_rP{pRow});
end

MPRow = M(row, pRows) * P(pRows, s_r);
 
end
