function slsOuts_alt = find_alt_impl_precise(sys, slsParams, slsOuts, Tc)
% Find alternate implementation, returned in slsOuts_alt
% Rearranged [Rc; Mc] = [R; M] [zI-A -B] [Rc; Mc] into the form
%   F2*[Rc; Mc] = F1 where F1, F2 do not depend on [Rc, Mc]
% Outputs
%    slsOuts_alt : contains clnorm and new R, M
% Inputs
%    sys         : LTISystem containing system matrices
%    slsParams   : SLSParams containing parameters
%    slsOuts     : contains info from SLS (original R, M)
%    Tc          : length of the approximate solution

statusTxt = 'Finding *exact* alt implementation';
statusTxt = [char(10), statusTxt, sprintf(' with Tc=%d', Tc)];
disp(statusTxt);

cvx_begin
cvx_solver sdpt3
cvx_precision low

variable Rcs(sys.Nx * (Tc-1), sys.Nx) % Rc{1} == 1 so it's not a variable
variable Mcs(sys.Nu * Tc, sys.Nx)

% helper function to convert cells of block matrices into giant matrix
% copied from get_F
get_range = @(idx, size) (size*(idx-1)+1:size*(idx-1)+size);

objective = 0;
for t = 1:Tc
    tx    = get_range(t-1, sys.Nx);
    tu    = get_range(t, sys.Nu);
    Mc{t} = Mcs(tu,:);

    if t==1
        Rc{t} = eye(sys.Nx);
    else        
        Rc{t} = Rcs(tx,:);
    end
    % L1 norm to enforce sparsity
    objective = objective + norm([sys.C1, sys.D12]*[Rc{t}; Mc{t}], 1);
end

% constraints
F  = get_F(sys, slsParams, slsOuts, Tc);
F1 = F(:, 1:sys.Nx);
F2 = F(:,sys.Nx+1:end);

F2 * [Rcs; Mcs] == -F1;

minimize(objective);
cvx_end

% outputs
slsOuts_alt.R_           = Rc;
slsOuts_alt.M_           = Mc;
slsOuts_alt.clnorm_      = objective;
slsOuts_alt.solveStatus_ = cvx_status;
