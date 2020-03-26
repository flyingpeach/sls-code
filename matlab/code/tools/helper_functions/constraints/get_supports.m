function [RSupp, MSupp, count] = get_supports(sys, params)
% Gets supports on R, M based on specified delay and/or localization

% default constraints (equivalent to having no constraints)
actDelay  = 0;
commSpeed = sys.Nx+1;
locality  = sys.Nx+1;
      
for i=1:length(params.constraints_)
    this_constr = params.constraints_{i};
    constrType  = this_constr(1);
    constrVal   = this_constr(2);
    
    switch constrType
        case SLSConstraint.ActDelay
            actDelay = constrVal;
        case SLSConstraint.CommSpeed
            commSpeed = constrVal;
        case SLSConstraint.Locality
            locality = constrVal;
    end
end

RSupp = cell(params.T_, 1); 
MSupp = cell(params.T_, 1);

commsAdj   = abs(sys.A) > 0;
localRSupp = commsAdj^(locality-1) > 0;

count = 0;
for t = 1:params.T_
    delayRSupp = commsAdj^(floor(max(0, commSpeed*(t - actDelay))));
    RSupp{t} = delayRSupp & localRSupp; % bitwise and   
    
    MSupp{t} = (abs(sys.B2)' * RSupp{t}) > 0; 
    count    = count + sum(sum(RSupp{t})) + sum(sum(MSupp{t}));
end

end