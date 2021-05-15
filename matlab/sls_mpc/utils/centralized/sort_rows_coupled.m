function [rCp, rUcp, nValsCp] = sort_rows_coupled(r, cpIdx)

len_r = length(r);

rCp     = cell(len_r, 1);
rUcp    = cell(len_r, 1);
nValsCp = 0; % number of rows with coupling

% Identify rows with coupling
for i = 1:len_r
    for j = 1:length(r{i})
        row = r{i}{j}; 
        if length(cpIdx{row}) <= 1 % there is only "self-coupling" (i.e. no coupling)
            rUcp{i}{end+1} = r{i}{j};
        else % there is coupling
            nValsCp = nValsCp + 1;            
            rCp{i}{end+1} = r{i}{j};
        end
    end
end

end