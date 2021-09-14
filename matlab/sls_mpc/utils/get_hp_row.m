function [HPRow, s_r] = get_hp_row(H, P, s_rP, rows, rows_neighbor, rowH)
% Returns row of H*Psi or H*Phi
% rows are rows of phi/psi that the current processor is in charge of
% P can be either Phi or Psi; s_rP can be s_rPhi or s_rPsi
% rows_cp are rows of psi that other processors are in charge of
%         that the current processor can access

pRows = find(H(rowH, :)); % all rows accessed via H*P
 
% sanity check that this processor has access to the desired row in Psi
if ~all(ismember(pRows, [rows, rows_neighbor]))
    mpc_error('Tried to access a row of Phi/Psi we do not have access to');
end

% Different rows will have different s_rPsi, take the union
s_r = [];

for pRow = pRows
    s_r = union(s_r, s_rP{pRow});
end

HPRow = H(rowH, pRows) * P(pRows, s_r);
 
end
