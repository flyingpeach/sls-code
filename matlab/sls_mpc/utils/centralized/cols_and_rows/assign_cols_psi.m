function c = assign_cols_psi(sys, T)
% c{i} represents the set of columns subsystem i solves for

Nx = sys.Nx;
c  = cell(Nx, 1);

for i = 1:Nx
    c{i} = zeros(1, T);
    for j = 1:T
        c{i}(j) = (j-1)*Nx + i;
    end
end

end