function M = build_from_cols(c, s_c, M_cols, M_size)

M = zeros(M_size);

for i = 1:length(c)
    for col = c{i}
        M(s_c{col}, col) = M_cols{col};
    end
end

end