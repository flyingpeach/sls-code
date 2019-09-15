function solnSpaceSize = get_soln_size(sys, slsParams, slsOuts, Tc, tol)
% Calculate size of solution space for each Tc
% Outputs
%     solnSpaceSize : the size of the solution space
% Inputs
%     sys           : LTISystem containing system matrices
%     slsParams     : SLSParams containing parameters
%     slsOuts       : contains info from SLS (original R, M)
%     Tc            : length of the approximate solution
%     tol           : tolerance to use for rank calculation

solnSpaceSize   = 0;
F      = get_F(sys, slsParams, slsOuts, Tc);
rankF  = rank(F, tol);
rankF2 = rank(F(:,sys.Nx+1:end), tol);

if rankF == rankF2
    solnSpaceSize = (Tc*(sys.Nx + sys.Nu - 1) - rankF) * sys.Nx;
end
    

