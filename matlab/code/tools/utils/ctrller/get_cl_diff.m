function [RDiff, MDiff] = get_cl_diff(sys, clMaps, ctrller, tTotal)
% Returns normalized difference between closed loop maps (R, M) and 
% closed loop maps (R_, M_) implemented by the ctrller (Rc, Mc)
% This will be inaccurate if unstable or tTotal too small

R = clMaps.R_; M = clMaps.M_;

T = length(R);

% get implemented cl maps
[R_, M_] = get_cl_map(sys, ctrller, tTotal);

RDiff = 0;
MDiff = 0;
RTotalNorm = 0;
MTotalNorm = 0;

for k=1:T
    RTotalNorm = RTotalNorm + norm(full(R{k}), 'fro');
    MTotalNorm = MTotalNorm + norm(full(M{k}), 'fro');
    RDiff = RDiff + norm(full(R{k} - R_{k}), 'fro');
    MDiff = MDiff + norm(full(M{k} - M_{k}), 'fro');
end
for k=T+1:tTotal % R, M =0
    RDiff = RDiff + norm(full(R_{k}), 'fro');
    MDiff = MDiff + norm(full(M_{k}), 'fro');
end

RDiff = RDiff / RTotalNorm;
MDiff = MDiff / MTotalNorm;


