function M_cols = separate_cols(c, s_c, M)

nCols  = size(M, 2);
M_cols = cell(nCols, 1);
                
for i = 1:length(c)
    for j = 1:length(c{i})
        col = c{i}{j};
        M_cols{col} = M(s_c{i}{j}, col);
    end
end

end