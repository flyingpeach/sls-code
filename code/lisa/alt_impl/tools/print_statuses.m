function statusTxt = print_statuses(sweepParams, slsOuts_alts)
% Print cvx statuses for each new CL implementation found
% Outputs
%     statusTxt    : contains list of swept parameters and corresponding status
% Inputs
%     sweepParams  : list of parameters (i.e. Tc, clDiffPen) corresponding
%                    to each slsOuts_alts{i}; should be same length as
%                    slsOuts_alts
%     slsOuts_alts : alternate CL implementations

statusTxt = [];


for i=1:length(sweepParams)
    param  = sweepParams(i);
    status = slsOuts_alts{i}.solveStatus_;

    statusTxt = [statusTxt, char(10), sprintf('%0.2e: %s', param, status)];
end
    

