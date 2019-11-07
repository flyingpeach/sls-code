function optCtrller = find_ctrller(sys, slsParams, slsOuts, Tc, settings)
% Find alternate implementation, returned in slsOuts_alt
%    optCtrller : contains (Rc, Mc) (controller implementation matrices)
% Inputs
%    sys        : LTISystem containing system matrices
%    slsParams  : SLSParams containing parameters
%    slsOuts    : contains info from SLS (original R, M)
%    Tc         : desired time horizon length of the approximate solution
%    settings   : AltImplSettings containing info on how to solve for alt impl

statusTxt = settings.sanity_check();
statusTxt = [char(10), sprintf('Finding alt implementation, %s, Tc=%d', statusTxt, Tc)];
disp(statusTxt);

F  = get_F(sys, slsParams, slsOuts, Tc);
F1 = F(:, 1:sys.Nx);
F2 = F(:,sys.Nx+1:end);

RMc_p         = F2 \ (-F1);
nullsp        = get_nullsp(F2, settings.eps_nullsp_);
solnSpaceSize = size(nullsp, 2);

cvx_begin
cvx_solver sdpt3
cvx_precision low

variable alpha(solnSpaceSize, sys.Nx)
RMc       = RMc_p + nullSp*alpha;
objective = norm(RMc, 1); % L1 norm minimization

[Rc, Mc] = block_to_cell(RMc, Tc, sys);

switch settings.mode_
    case OptL1Delayed
    % todo
    



%       OptL1 % Minimizes L1 norm of [Rc; Mc] 
%       OptL1Delayed   % OptL1 with delay constraints only
%       OptL1Localized % OptL1 with locality constraints only
%       OptL1DAndL     % OptL1 with delay and locality constraints
%       EncourageDelay % OptL1 that encourages tolerance for communication delay
%       EncourageLocal % OptL1 that encourages tolerance for locality

end


% local functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function range = get_range(idx, size)
% helper function to convert cells of block matrices into giant matrix
% copied from get_F (inline)
range = size*(idx-1)+1:size*(idx-1)+size;
end 


function [Rc, Mc] = block_to_cell(RMc, Tc, sys)
% Rc, Mc are cell structure (i.e. Rc{t})
% RMc is one stacked matrix (Rc first, then Mc)
Rcs = RMc(1:sys.Nx * (Tc-1), :);
Mcs = RMc(sys.Nx * (Tc-1) + 1:end, :);

for t = 1:Tc
    tx    = get_range(t-1, sys.Nx);
    tu    = get_range(t, sys.Nu);
    Mc{t} = Mcs(tu,:);

    if t==1
        Rc{t} = eye(sys.Nx);
    else        
        Rc{t} = Rcs(tx,:);
    end
end
end
