function ConstrMtx = build_constr_mtx(sys, params)

Nx   = sys.Nx;
tFIR = params.tFIR_;
K    = params.constrMtx_;

ConstrMtx = [zeros(size(K)) zeros(2*Nx,(tFIR-1)*Nx)];
ConstrMtx = [ConstrMtx; zeros(2*Nx,Nx) zeros(size(K)) zeros(2*Nx,(tFIR-2)*Nx)];

for t = 2:tFIR-1
    ConstrMtx = [ConstrMtx; zeros(2*Nx,t*Nx) K zeros(2*Nx,(tFIR-t-1)*Nx)];
end

end