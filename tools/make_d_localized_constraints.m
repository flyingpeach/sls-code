function [Rsupport, Msupport, count] = make_d_localized_constraints(sys, params)
% Inputs
%    sys    : an LTISystem
%    params : SLSParams containing parameters

commsAdj  = abs(sys.A) > 0;
localityR = commsAdj^(params.d_-1) > 0;

count = 0;
for t = 1:params.tFIR_
    Rsupport{t} = min(commsAdj^(floor(max(0, params.cSpeed_*(t-params.actDelay_)))),localityR) > 0;
    Msupport{t} = (abs(sys.B2)'*Rsupport{t}) > 0;
    count       = count + sum(sum(Rsupport{t})) + sum(sum(Msupport{t}));
end