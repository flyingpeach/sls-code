function statusTxt = print_statuses(slsOuts_alts)
% Print cvx statuses for each new CL implementation found
% Outputs
%     statusTxt    : contains list of Tc and corresponding status
% Inputs
%     slsOuts_alts : alternate CL implementations (one per different Tc)

statusTxt = [];

for i=1:length(slsOuts_alts)
    Tc = length(slsOuts_alts{i}.R_);

    status    = slsOuts_alts{i}.solveStatus_;
    statusTxt = [statusTxt, char(10), sprintf('Tc=%d, %s', Tc, status)];     
end
    

