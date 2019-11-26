function ctrller = find_ctrller(sys, slsOuts, cParams)
% Find alternate implementation, returned in slsOuts_alt
%    ctrller   : Ctrller containing (Rc, Mc)
% Inputs
%    sys       : LTISystem containing system matrices
%    slsOuts   : contains info from SLS (original R, M)
%    cParams   : CtrllerParams containing info on how to solve for alt impl

statusTxt = cParams.sanity_check();
statusTxt = [char(10), sprintf('Finding a controller, %s', statusTxt)];
disp(statusTxt);

F  = get_ctrller_constraint(sys, slsOuts, cParams.T_);
F1 = F(:, 1:sys.Nx);
F2 = F(:,sys.Nx+1:end);

RMc_p         = F2 \ (-F1);
nullsp        = get_nullsp(F2, cParams.eps_nullsp_);
solnSpaceSize = size(nullsp, 2);

cvx_begin
cvx_solver sdpt3
cvx_precision low

variable nullCoeff(solnSpaceSize, sys.Nx)
RMc_free = RMc_p + nullsp*nullCoeff; % Rc and Mc without constraints

% modes that have no constraints (other than basic Rc/Mc constr)
freeModes   = [SLSMode.Basic, SLSMode.EncourageDelay, SLSMode.EncourageLocal];

% modes that have locality / delay constraints
constrModes = [SLSMode.Delayed, SLSMode.Localized, SLSMode.DAndL];

if ismember(cParams.mode_, freeModes)
    objective = norm(RMc_free, 1); % L1 norm minimization
    
    [Rc, Mc] = block_to_cell(RMc_free, cParams, sys);
    if cParams.mode_ == SLSMode.EncourageDelay
        for t=1:cParams.T_
            % formulating these constraints are slow, so output progress
            disp(sprintf('Adding objectives for %d of %d matrices', t, cParams.T_));
        
            BM = sys.B2 * Mc{t};
            for i=1:sys.Nx
                for j=1:sys.Nx % TODO: this distance metric only for chain
                    objective = objective + cParams.fastCommPen_ * exp(abs(i-j)-t)*norm(Rc{t}(i,j));
                    objective = objective + cParams.fastCommPen_ * exp(abs(i-j)-t)*norm(BM(i,j));
                end
            end
        end
    elseif cParams.mode_ == SLSMode.EncourageLocal
        for t=1:cParams.T_
            % formulating these constraints are slow, so output progress
            disp(sprintf('Adding objectives for %d of %d matrices', t, cParams.T_));

            BM = sys.B2 * Mc{t};
            for i=1:sys.Nx
                for j=1:sys.Nx % TODO: this distance metric only for chain
                    objective = objective + cParams.nonLocalPen_ * exp(abs(i-j))*norm(Rc{t}(i,j));
                    objective = objective + cParams.nonLocalPen_ * exp(abs(i-j))*norm(BM(i,j));
                end
            end
        end
    end
elseif ismember(cParams.mode_, constrModes)
    [RSupp, MSupp, count] = get_supports(sys, cParams);  
    nRowsRMc = sys.Nx*(cParams.T_-1) + sys.Nu*cParams.T_;
    expression RMc_constr(nRowsRMc, sys.Nx)
    variable RMSupp(count)

    spot = 0;
    for t=1:cParams.T_
        if t > 1
            [is, js] = find(RSupp{t}); % rows, cols of support
            num      = length(is);     % # nonzero elements
            is       = sys.Nx * (t-2) + is; % shift row
            RMc_constr(is + (js-1)*nRowsRMc) = RMSupp(spot+1:spot+num);
            spot = spot + num;
        end
        [is, js] = find(MSupp{t});
        num      = length(is);
        is       = sys.Nx * (cParams.T_-1) + sys.Nu * (t-1) + is;
        RMc_constr(is + (js-1)*nRowsRMc) = RMSupp(spot+1:spot+num);
        spot = spot + num;
    end
    
    % heuristic weighted distance from CL map plus L1 objective
    objective = cParams.CLDiffPen_*norm(RMc_free - RMc_constr, 2) + norm(RMc_constr, 1); 

    [Rc, Mc] = block_to_cell(RMc_constr, cParams, sys);
end

minimize(objective);
cvx_end

ctrller = Ctrller(); % output
ctrller.Rc_ = Rc;
ctrller.Mc_ = Mc;
end


% local functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function range = get_range(idx, size)
% helper function to convert cells of block matrices into giant matrix
% copied from get_F (inline)
range = size*(idx-1)+1:size*(idx-1)+size;
end 


function [Rc, Mc] = block_to_cell(RMc, cParams, sys)
% Rc, Mc are cell structure (i.e. Rc{t})
% RMc is one stacked matrix (Rc first, then Mc)
Rcs = RMc(1:sys.Nx * (cParams.T_-1), :);
Mcs = RMc(sys.Nx * (cParams.T_-1) + 1:end, :);

for t = 1:cParams.T_
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
