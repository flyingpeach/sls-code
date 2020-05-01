function [Psi_rows, Lambda_rows] = separate_rows_grouped(sys, r, s_r, Psi, Lambda)

Nx = sys.Nx;

Psi_rows    = cell(1, Nx);
Lambda_rows = cell(1, Nx);

for i = 1:Nx
    % a bit hacky: s_r{i}(j) is equal for all j, so just take j=1
    sri = s_r{i}(1, :); 
    Psi_rows{i}    = Psi(r{i}, sri);
    Lambda_rows{i} = Lambda(r{i}, sri);
end

end