function XiSupp = get_sparsity_xi(H, G, Nx, PsiSupp)
% PsiSupp is optional:
% if unspecified, we will use H = Xi*G to calculate sparsity
% if specified,   we will use H*Psi{2:T} = Xi*G to calculate sparsity

numRowsH = size(H, 1);
numRowsG = size(G, 1);
numColsG = size(G, 2); % equal to Nx*(T-1)
XiSupp   = true(numRowsH, numRowsG); % start with full support

Jis = cell(numRowsH, 1); % for row i of H, the set of j s.t. H(i,j) ~= 0
Lks = cell(numColsG, 1); % for col k of G, the set of l s.t. G(l,k) ~= 0

for i=1:numRowsH
    Jis{i} = find(H(i,:));
end

for k=1:numColsG
    Lks{k} = find(G(:,k));
end

mtx = eye(Nx);
if nargin == 4
    mtx = PsiSupp(:, Nx+1:end);
end

for i=1:numRowsH
    for k=1:numColsG       
        if all(mtx(Jis{i}, k) == 0)
            XiSupp(i, Lks{k}) = false;
        end
    end
end

end
