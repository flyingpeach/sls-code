function ConstrMtx = build_constr_mtx(sys, params)

Nx   = sys.Nx;
tFIR = params.tFIR_;

% TODO: hacky; ignore input constraints (use zeros)
K = params.stateConsMtx_;

K_ = zeros(2*Nx,Nx); j = 0;
for i = 1:2*Nx
    if mod(i,4) == 1 || mod(i,4) == 2
        j = j + 1;
        K_(i,:) = K(j,:);
    else
    end
end

ConstrMtx = [zeros(size(K_)) zeros(2*Nx,(tFIR-1)*Nx)];
ConstrMtx = [ConstrMtx; zeros(2*Nx,Nx) zeros(size(K_)) zeros(2*Nx,(tFIR-2)*Nx)];
for t = 2:tFIR-1
    ConstrMtx = [ConstrMtx; zeros(2*Nx,t*Nx) K_ zeros(2*Nx,(tFIR-t-1)*Nx)];
end

end