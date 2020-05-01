function M_cols = separate_cols(sys, c, s_c, M)

Nx     = sys.Nx;
M_cols = cell(Nx, 1);

for i = 1:Nx
    M_cols{i} = M(s_c{i}, c{i});
end
        
