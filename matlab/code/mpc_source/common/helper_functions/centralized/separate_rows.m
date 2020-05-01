function [Psi_rows, Lambda_rows] = separate_rows(sys, tFIR, r, s_r, Psi, Lambda)

Nx = sys.Nx; Nu = sys.Nu;
nRows = Nx*tFIR + Nu*(tFIR-1);

Psi_rows    = cell(nRows, 1);
Lambda_rows = cell(nRows, 1);
                
for i = 1:Nx
    for j = 1:length(r{i})
        row = r{i}(j);    
        Psi_rows{row}    = Psi(row, s_r{i}(j, :));
        Lambda_rows{row} = Lambda(row, s_r{i}(j, :));
    end
end

end