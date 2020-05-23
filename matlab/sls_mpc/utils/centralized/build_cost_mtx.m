function CostMtx = build_cost_mtx(params)

tFIR = params.tFIR_;

CostMtx = [];
for t = 0:tFIR-1
    CostMtx = blkdiag(CostMtx, params.QSqrt_);
end
for t = 0:tFIR-2
    CostMtx = blkdiag(CostMtx, params.RSqrt_);
end    

end