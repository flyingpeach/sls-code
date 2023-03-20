function [mtx, parTime] = get_local_subspace(sys, x0, params, eps)
% We use tFIR and locality size from params
% Space of local trajectories is proportional to the size of Image(mtx)
% eps : how to determine whether solution exists

Nx   = sys.Nx; Nu = sys.Nu; T = params.tFIR_;
nPhi = Nx*T + Nu*(T-1);
ZAB  = sparse(get_constraint_zab(sys, T));
IO   = [eye(Nx); zeros(Nx*(T-1), Nx)];
k    = vec(IO);

PsiSupp  = get_sparsity_psi(sys, params);
PsiSupp  = PsiSupp(:, 1:Nx); % The rest of the matrix is for robust only
suppIdx  = find(PsiSupp);

mtx = [];

parTimes = zeros(Nx, 1);
for i=1:Nx
    tic;
    idxStart = (i-1)*nPhi+1;
    idxEnd   = i*nPhi;
    myIdx    = suppIdx(suppIdx >= idxStart & suppIdx <= idxEnd) - (i-1)*nPhi;
    nCols    = length(myIdx);
    
    % Check if there is a solution
    Hi       = ZAB(:,myIdx);
    ki       = k((i-1)*Nx*T+1:i*Nx*T);    
    testSol  = Hi\ki;
    
    % Note: using max instead of 2-norm to make this check less
    %       dimension-dependent
    if max(abs(Hi*testSol - ki)) > eps % No solution exists
        mtx = 0;
        parTimes(i) = toc;
        parTime = mean(parTimes);
        return;
    end
    
    IHHi     = speye(nCols) - Hi\Hi;

    % Equivalent to above, but slower
    % Uncomment this out for numerical_example.m
    % IHHi     = speye(nCols) - pinv(full(Hi))*Hi; 

    if(any(isnan(IHHi(:)))) % Conditioning poor, no solution available
        mtx = 0;
        parTimes(i) = toc;
        parTime = mean(parTimes);
        return;
    end
    
    zeroCols = false(nCols,1); % Get rid of zero columns
    for j=1:nCols
        if isempty(find(IHHi(:,j),1))
            zeroCols(j) = true;
        end
    end
    
    % Comment this out for numerical_example.m (for full matrix)
    IHHi(:,zeroCols) = []; 
    
    xi  = x0(i)*speye(nPhi); 
    mtx = [mtx xi(Nx+1:end,myIdx)*IHHi];
    parTimes(i) = toc;
end

parTime = mean(parTimes);