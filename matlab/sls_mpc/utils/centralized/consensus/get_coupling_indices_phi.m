function cpIdx = get_coupling_indices_phi(C, K)
% cpIdx : per row of Phi, which other rows of phi are coupled to it
% C     : cost matrix (augmented to size of Phi)
% K     : constraint matrix (augmented to size of Phi)

nPhi  = length(C);
cpIdx = cell(1, nPhi); 

for i = 1:nPhi
    for j = 1:nPhi
        if (C(i,j) ~= 0 || K(i,j) ~=0)
            cpIdx{i} = unique([cpIdx{i} j]);
            cpIdx{j} = unique([cpIdx{j} i]);
        end
    end
end

end