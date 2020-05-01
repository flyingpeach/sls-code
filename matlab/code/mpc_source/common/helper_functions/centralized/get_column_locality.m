function [c, s_c] = get_column_locality(sys, tFIR, r_loc, m_loc)
% c{i}   represents the set of columns controller i solves for
% s_c{i} represents the set of rows associated to the columns in c{i} 

Nx = sys.Nx; 

rm_loc = []; % same size as Phi, Psi, Lambda
for t = 1:tFIR
    rm_loc = [rm_loc; r_loc];
end
for t = 1:tFIR-1
    rm_loc = [rm_loc; m_loc];
end

c     = cell(Nx, 1);
s_c   = cell(Nx, 1);

for i = 1:Nx
    c{i} = i; % each controller receives one column
    s_c{i} = find(rm_loc(:,i)); % find nonzero rows of this column    
end

end