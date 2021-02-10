function ConstrMtx = build_constr_mtx(sys, params)

tFIR = params.tFIR_;

ConstrMtx = [];

K1 = zeros(sys.Nx);
K2 = zeros(sys.Nu);
if params.has_state_cons()
    K1 = params.stateConsMtx_;
end
if params.has_input_cons()
    K2 = params.inputConsMtx_;
end
    
ConstrMtx = zeros(sys.Nx);
for t = 2:tFIR
    ConstrMtx = blkdiag(ConstrMtx, K1);
end
for t = 1:tFIR-1
    ConstrMtx = blkdiag(ConstrMtx, K2);
end    

end