function [RSupp, MSupp, count] = get_supports(sys, params)
% Gets supports on R, M based on specified delay and/or localization
% Expects that it's only called if mode involves delay and/or localization

if params.mode_ == SLSMode.Delayed % i.e. no locality constraint
    d = sys.Nx+1; % "locality" is entire graph; used to simplify code
                  % this is the same as having no locality constraint
else
    d = params.d_;
end

commsAdj = abs(sys.A) > 0;
localRSupp = commsAdj^(d-1) > 0;

count = 0;
for t = 1:params.tFIR_
    if params.mode_ == SLSMode.Localized % locality constraint only
        RSupp{t} = localRSupp > 0;
    else % include delay constraint as well
        delayRSupp = commsAdj^(floor(max(0, params.cSpeed_*(t-params.actDelay_))));
        RSupp{t} = min(delayRSupp, localRSupp) > 0;    
    end
    
    MSupp{t} = (abs(sys.B2)' * RSupp{t}) > 0; 
    count    = count + sum(sum(RSupp{t})) + sum(sum(MSupp{t}));    
end

end