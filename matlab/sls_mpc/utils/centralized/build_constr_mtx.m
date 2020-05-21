function ConstrMtx = build_constr_mtx(sys, params)

tFIR = params.tFIR_;

ConstrMtx = [];

% TODO: currently assumes K1 is Nx by Nx and somewhat diagonal
K1 = params.stateConsMtx_;

% TODO: currently ignores input constraint matrix
K2 = zeros(sys.Nu);   
    
for t = 0:tFIR-1
    ConstrMtx = blkdiag(ConstrMtx, K1);
end
for t = 0:tFIR-2
    ConstrMtx = blkdiag(ConstrMtx, K2);
end    

end