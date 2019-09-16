function slsOuts_alt = find_alt_impl_ls(sys, slsParams, slsOuts, Tc)
% Find alternate implementation via least squares; should be the alternate
% implementation with closest CL map to original R,M
% Note that this is meant to be used for overconstrained systems; an exact
% solution will be provided for a non-overconstrained system
% Outputs
%    slsOuts_alt : contains clnorm (currently L1 norm) and new R, M
% Inputs
%    sys         : LTISystem containing system matrices
%    slsParams   : SLSParams containing parameters
%    slsOuts     : contains info from SLS (original R, M)
%    Tc          : length of the approximate solution

statusTxt = 'Finding least squares alt implementation';
statusTxt = [char(10), statusTxt, sprintf(' with Tc=%d', Tc)];
disp(statusTxt);

F  = get_F(sys, slsParams, slsOuts, Tc);
F1 = F(:, 1:sys.Nx);
F2 = F(:,sys.Nx+1:end);

% F2 * [Rcs; Mcs] == -F1;
RMc = F2 \ (-F1); % if overconstrained, least squares solution

% convert to cell form for easier calculation
Rcs = RMc(1:sys.Nx * (Tc-1), :);
Mcs = RMc(sys.Nx * (Tc-1) + 1:end, :);

% helper function to convert cells of block matrices into giant matrix
% copied from get_F
get_range = @(idx, size) (size*(idx-1)+1:size*(idx-1)+size);

L1norm = 0;
for t = 1:Tc
    tx    = get_range(t-1, sys.Nx);
    tu    = get_range(t, sys.Nu);
    Mc{t} = Mcs(tu,:);

    if t==1
        Rc{t} = eye(sys.Nx);
    else        
        Rc{t} = Rcs(tx,:);
    end
    L1norm = L1norm + norm([sys.C1, sys.D12]*[Rc{t}; Mc{t}], 1);
end

slsOuts_alt.R_           = Rc;
slsOuts_alt.M_           = Mc;
slsOuts_alt.clnorm_      = L1norm;
slsOuts_alt.solveStatus_ = 'Analytic';
end