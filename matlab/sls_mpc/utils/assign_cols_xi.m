function c = assign_cols_xi(sys, T)
% c{i} represents the set of columns subsystem i solves for

% Without coupling, we will have a UB and LB on each disturbance
% so G will have 2*Nx*(T-1) rows; Xi will have this many columns

Nx = sys.Nx;
c  = cell(Nx, 1);

for i = 1:Nx
    c{i} = zeros(1, 2*(T-1));
    for j = 1:2*(T-1)
        c{i}(j) = Nx*(j-1) + i;
    end
end

end