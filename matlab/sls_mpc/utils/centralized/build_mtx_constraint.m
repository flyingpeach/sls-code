function ConstrMtx = build_mtx_constraint(sys, params)

T = params.tFIR_;

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
for t = 2:T
    ConstrMtx = blkdiag(ConstrMtx, K1);
end
for t = 1:T-1
    ConstrMtx = blkdiag(ConstrMtx, K2);
end    

end