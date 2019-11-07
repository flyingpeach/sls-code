function ctrller = find_ctrller(sys, slsParams, slsOuts, cParams)
% Find alternate implementation, returned in slsOuts_alt
%    ctrller   : Ctrller containing (Rc, Mc)
% Inputs
%    sys       : LTISystem containing system matrices
%    slsParams : SLSParams containing parameters
%    slsOuts   : contains info from SLS (original R, M)
%    cParams   : CtrllerParams containing info on how to solve for alt impl

statusTxt = cParams.sanity_check();
statusTxt = [char(10), sprintf('Finding a controller, %s', statusTxt)];
disp(statusTxt);

F  = get_F(sys, slsParams, slsOuts, cParams.tc_);
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
        for t=1:cParams.tc_
            % formulating these constraints are slow, so output progress
            [char(10), sprintf('Adding objectives for %d of %d matrices', t, cParams.tc_)]
        
            BM = sys.B2 * Mc{t};
            for i=1:sys.Nx
                for j=1:sys.Nx % note: this distance metric only for chain
                    objective = objective + cParams.fascParams.tc_ommPen_ * exp(abs(i-j)-t)*norm(Rc{t}(i,j));
                    objective = objective + cParams.fascParams.tc_ommPen_ * exp(abs(i-j)-t)*norm(BM(i,j));
                end
            end
        end
    elseif cParams.mode_ == SLSMode.EncourageLocal
        for t=1:cParams.tc_
            % formulating these constraints are slow, so output progress
            [char(10), sprintf('Adding objectives for %d of %d matrices', t, cParams.tc_)]
        
            BM = sys.B2 * Mc{t};
            for i=1:sys.Nx
                for j=1:sys.Nx % note: this distance metric only for chain
                    objective = objective + cParams.nonLocalPen_ * exp(abs(i-j))*norm(Rc{t}(i,j));
                    objective = objective + cParams.nonLocalPen_ * exp(abs(i-j))*norm(BM(i,j));
                end
            end
        end
    end
elseif ismember(cParams.mode_, constrModes)
    expression RMc_constr((cParams.tc_-1)*sys.Nx + cParams.tc_*sys.Nu, sys.Nx)

    [Rc, Mc] = block_to_cell(RMc_constr, cParams, sys);
    [RSupp, MSupp, count] = get_supports(sys, cParams);  

    % TODO: copied from state_fdbk_sls.m
    variable RMSupp(count)
    spot = 0;
    for t = 1:params.tFIR_
        suppR = find(RSupp{t});
        num = sum(sum(RSupp{t}));
        Rc{t}(suppR) = RMSupp(spot+1:spot+num);
        spot = spot + num;

        suppM = find(MSupp{t});
        num = sum(sum(MSupp{t}));
        Mc{t}(suppM) = RMSupp(spot+1:spot+num);
        spot = spot + num;
    end

    % heuristic distance from CL map plus L1 objective
    objective = norm(RMc_free - RMc_constr, 2) + norm(RMc_constr, 1); 
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
Rcs = RMc(1:sys.Nx * (cParams.tc_-1), :);
Mcs = RMc(sys.Nx * (cParams.tc_-1) + 1:end, :);

for t = 1:cParams.tc_
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
