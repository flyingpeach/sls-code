function statusTxt = print_statuses(sys, slsParams, slsOuts, slsOuts_alts, tol)
% Print cvx statuses for each new CL implementation found
% Outputs
%     statusTxt    : contains list of Tc and corresponding status
% Inputs
%     sys          : LTISystem containing system matrices
%     slsParams    : SLSParams containing parameters
%     slsOuts      : contains info from SLS (original R, M)
%     slsOuts_alts : alternate CL implementations (one per different Tc)
%     tol          : tolerance to use for rank calculation

statusTxt = [];

for i=1:length(slsOuts_alts)
    Tc = length(slsOuts_alts{i}.R_);

    F  = get_F(sys, slsParams, slsOuts, Tc);
    F2 = F(:,sys.Nx+1:end);
    
    solnSpaceSize = size(get_soln_sp(F2, tol), 2);
    status        = slsOuts_alts{i}.solveStatus_;

    statusTxt = [statusTxt, char(10), sprintf('Tc=%d, %s, solnSpaceSize=%d', Tc, status, solnSpaceSize)];     
end
    

