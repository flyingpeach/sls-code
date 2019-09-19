function [RSupp, MSupp, count] = get_localized_supports(sys, params)
% TODO: copied from state_fdbk_sls
commsAdj  = abs(sys.A) > 0;
localityR = commsAdj^(params.d_-1) > 0;

count = 0;
for t = 1:params.tFIR_
    RSupp{t} = min(commsAdj^(floor(max(0, params.cSpeed_*(t-params.actDelay_)))),localityR) > 0;
    MSupp{t} = (abs(sys.B2)'*RSupp{t}) > 0;
    count       = count + sum(sum(RSupp{t})) + sum(sum(MSupp{t}));
end
end
