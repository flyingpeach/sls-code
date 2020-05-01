function M = build_from_rows(sys, r, s_r, M_rows, M_size)

M = zeros(M_size);

for i = 1:sys.Nx
    for j = 1:length(r{i})
        row = r{i}(j);
        M(row, s_r{i}{j}) = M_rows{row};
    end
end

end