function M_rows = separate_rows(r, s_r, M)

nRows  = size(M, 1);
M_rows = cell(nRows, 1);
                
for i = 1:length(r)
    for j = 1:length(r{i})
        row = r{i}{j};
        M_rows{row} = M(row, s_r{i}{j});
    end
end

end