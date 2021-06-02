function c = assign_cols_phi(sys)
% c{i}   represents the set of columns subsystem i solves for

Nx = sys.Nx;

c = cell(sys.Nx, 1);
for i = 1:Nx
    c{i} = [i]; % each controller receives one column
end

end