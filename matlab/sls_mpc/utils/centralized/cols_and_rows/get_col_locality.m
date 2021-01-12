function s_c = get_col_locality(c, SuppMtx)
% c{i}   represents the set of columns controller i solves for
% s_c{i} represents the set of rows associated to the columns in c{i}

Nx  = length(c);
s_c = cell(Nx, 1);

for i=1:Nx
    numColsi = length(c{i});
    s_c{i}   = cell(numColsi, 1);

    for j = 1:numColsi
        col       = c{i}{j};
        s_c{i}{j} = find(SuppMtx(:,col));
    end
end

end