function [c, s_c] = get_col_locality(sys, PhiSupp)
% c{i}   represents the set of columns subsystem i solves for
% s_c{i} represents the set of rows associated to the columns in c{i} 

for i = 1:sys.Nx
    c{i} = i; % each controller receives one column
    s_c{i} = find(PhiSupp(:,i)); % find nonzero rows of this column    
end

end