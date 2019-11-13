function solnSpaceSize = check_soln_space_size(sys, slsOuts, Tc, tol)
% Returns the size of the solution space to F2[Rc; Mc] = -F1
% If no solution, return 0
% Outputs
%    solnSpaceSize : size of the solution space
% Inputs
%    sys           : LTISystem containing system matrices
%    slsParams     : SLSParams containing parameters
%    slsOuts       : contains info from SLS (original R, M)
%    Tc            : length of the approximate solution
%    tol           : tolerance to use for rank / nullspace calculations

F  = get_F(sys, slsOuts, Tc);
F2 = F(:,sys.Nx+1:end);

solnSpaceSize = 0;
if rank(F, tol) == rank(F2, tol)
    nullsp        = get_nullsp(F2, tol);
    solnSpaceSize = size(nullsp, 2);
end

