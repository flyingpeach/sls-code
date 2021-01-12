function s_r = get_row_locality(r, SuppMtx)
% r{i}   represents the set of rows controller i solves for
% s_r{i} represents the set of columns associated to the rows in r{i} 

Nx  = length(r);
s_r = cell(Nx, 1);

for i=1:Nx
    numRowsi = length(r{i});
    s_r{i}   = cell(numRowsi, 1);

    for j = 1:numRowsi
        row       = r{i}{j};
        s_r{i}{j} = find(SuppMtx(row,:));
    end
end

end