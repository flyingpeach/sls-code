function [CPsiRow, s_r] = get_cpsi_row(C, Psi, s_rPsi, rows, rows_neighbor, row)
% Returns row of C*Psi
% rows are rows of phi/psi that the current processor is in charge of
% s_rPsi is the sparsity of Psi we want to use
% rows_cp are rows of psi that other processors are in charge of
%         that the current processor can access

psiRows = find(C(row, :)); % all rows accessed via C*Psi
 
% sanity check that this processor has access to the desired row in Psi
if ~all(ismember(psiRows, [rows, rows_neighbor]))
    mpc_error('Tried to access a row of Psi we do not have access to');
end

% Different rows will have different s_rPsi, take the union
s_r = [];

for psiRow = psiRows
    s_r = union(s_r, s_rPsi{psiRow});
end

CPsiRow = C(row, psiRows) * Psi(psiRows, s_r);
 
end
