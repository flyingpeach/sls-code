function [rCp, rUcp, row2cp, nCp] = sort_coupled_rows(r, cpIdx)
% rCp     : per subsystem, rows that are coupled 
% rUcp    : per subsystem, rows that are not coupled
% nValsCp : total number of coupled values
% row2cp(i) takes the ith row of the Phi matrix and maps it to
%     the index in the coupled variables (X, Y, Z). If the ith row is not
%     coupled, row2cp(i) will give 0

nSubsys = length(r);
rCp     = cell(nSubsys, 1);
rUcp    = cell(nSubsys, 1);
cp2row  = [];

% Identify rows with coupling
for i = 1:nSubsys
    for row = r{i}
        if length(cpIdx{row}) <= 1 % no coupling
            rUcp{i}(end+1) = row;
        else % there is coupling
            rCp{i}(end+1) = row;
            cp2row(end+1) = row;
        end
    end
end

nCp = length(cp2row);
row2cp  = zeros(1, nCp);
for i = 1:nCp
    row2cp(cp2row(i)) = i;
end

end