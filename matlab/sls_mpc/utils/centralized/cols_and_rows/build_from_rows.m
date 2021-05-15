function M = build_from_rows(r, s_r, M_rows, M_size)

M = zeros(M_size);

for i = 1:length(r)
    for row = r{i}
        M(row, s_r{row}) = M_rows{row};
    end
end

end