function statusTxt = print_statuses(sys, slsOuts_alts, rankFs, rankF2s)
% Print cvx statuses for each new CL implementation found
% Outputs
%     statusTxt    : contains list of Tc and corresponding status
% Inputs
%     sys          : LTISystem containing system matrices
%     slsOuts_alts : alternate CL implementations (one per different Tc)
%     rankFs       : ranks of F (one per different Tc)
%     rankF2s      : ranks of F2 (one per different Tc)

statusTxt = [];

for i=1:length(slsOuts_alts)
    Tc = length(slsOuts_alts{i}.R_);

    solnSpaceSize = 0;
    if rankFs(i) == rankF2s(i);
        solnSpaceSize = (Tc*(sys.Nx + sys.Nu - 1) - rankFs(i)) * sys.Nx;
    end

    status    = slsOuts_alts{i}.solveStatus_;
    statusTxt = [statusTxt, char(10), sprintf('Tc=%d, %s, solnSpaceSize=%d', Tc, status, solnSpaceSize)];     
end
    

