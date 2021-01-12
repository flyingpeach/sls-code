function M_rows = separate_rows(sys, tFIR, r, s_r, M)

Nx = sys.Nx; Nu = sys.Nu;
nRows = Nx*tFIR + Nu*(tFIR-1);

M_rows = cell(nRows, 1);
                
for i = 1:Nx
    for j = 1:length(r{i})
        row = r{i}{j};
        M_rows{row} = M(row, s_r{i}{j});
    end
end

end