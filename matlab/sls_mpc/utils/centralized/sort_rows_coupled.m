function [rCp, rUcp, nValsCp] = sort_rows_coupled(r, cpIdx)

len_r = length(r);

rCp     = cell(len_r, 1);
rUcp    = cell(len_r, 1);
nValsCp = 0; % number of rows with coupling

% Identify rows with coupling
for i = 1:len_r
    for row = r{i}
        if length(cpIdx{row}) <= 1 % no coupling
            rUcp{i}(end+1) = row;
        else % there is coupling
            nValsCp = nValsCp + 1;            
            rCp{i}(end+1) = row;
        end
    end
end

end