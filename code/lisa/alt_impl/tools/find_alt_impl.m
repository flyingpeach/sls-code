function slsOuts_alt = find_alt_impl(sys, slsParams, slsOuts, Tc, settings)
% Find alternate implementation, returned in slsOuts_alt
%    slsOuts_alt : contains clnorm and Rc, Mc
% Inputs
%    sys         : LTISystem containing system matrices
%    slsParams   : SLSParams containing parameters
%    slsOuts     : contains info from SLS (original R, M)
%    Tc          : length of the approximate solution
%    settings  : AltImplSettings containing info on how to solve for alt impl

% TODO: this output is a misnomer as it's not really an "slsout"

statusTxt = settings.sanity_check();
statusTxt = [char(10), sprintf('Finding alt implementation, %s, Tc=%d', statusTxt, Tc)];
disp(statusTxt);

F  = get_F(sys, slsParams, slsOuts, Tc);
F1 = F(:, 1:sys.Nx);
F2 = F(:,sys.Nx+1:end);

if settings.mode_ == AltImplMode.Analytic
    slsOuts_alt = fai_analytic(sys, Tc, F1, F2);   
elseif settings.mode_ ~= AltImplMode.ApproxLS
    tol = settings.tol_;
    if settings.mode_ ~= AltImplMode.ApproxLeaky
        if rank(F2, tol) ~= rank([F1 F2], tol)
            fprintf('Solution infeasible!');
            slsOuts_alt = 0;
            return
        end
    end
    slsOuts_alt = fai_nullsopt(sys, Tc, F1, F2, tol, settings);
else
    tol = settings.svThresh_;
    slsOuts_alt = fai_nullsopt(sys, Tc, F1, F2, tol, settings);
end
end


% local functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function range = get_range(idx, size)
% helper function to convert cells of block matrices into giant matrix
% copied from get_F (inline)
range = size*(idx-1)+1:size*(idx-1)+size;
end 


function [Rc, Mc] = block_to_cell(Rcs, Mcs, Tc, sys)
% Rc, Mc are cell structure (i.e. Rc{t}); Rcs, Mcs are stacked matrices
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


function slsOuts_alt = fai_nullsopt(sys, Tc, F1, F2, tol, settings)

RMc_p         = F2 \ (-F1); % particular solution
nullSp        = get_nullsp(F2, tol);
solnSpaceSize = size(nullSp, 2);

cvx_begin
cvx_solver sdpt3
cvx_precision low

variable coeff(solnSpaceSize, sys.Nx)
RMc = RMc_p + nullSp*coeff;
objective  = norm(RMc, 1);
    
if settings.mode_ == AltImplMode.ApproxLeaky
    variable RMcSlack(size(RMc_p, 1), sys.Nx)
    RMc = RMc + RMcSlack;
    objective = objective + settings.clDiffPen_ * norm(RMcSlack);
end

Rcs = RMc(1:sys.Nx * (Tc-1), :);
Mcs = RMc(sys.Nx * (Tc-1) + 1:end, :);

[Rc, Mc] = block_to_cell(Rcs, Mcs, Tc, sys);

if settings.mode_ == AltImplMode.StrictDelay
    [RZeros, MZeros] = get_delay_constraints(sys, Tc, settings.delay_);
    for t=1:Tc
        Rc{t}(RZeros{t}) == 0;
        if t > 1
            Mc{t}(MZeros{t}) == 0; % special case; let it leak, then zero
                                   % otherwise commonly infeasible
        end
    end
    objective = objective + settings.m1NonzeroPen_ * norm(Mc{1}(MZeros{1}));
end

if settings.mode_ == AltImplMode.StrictLocal
    [RZeros, MZeros] = get_locality_constraints(sys, settings.locality_);
    for t=1:Tc
        Rc{t}(RZeros) == 0;
        if t > 1
            Mc{t}(MZeros) == 0; % special case; let it leak
        end
    end
    objective = objective + settings.m1NonzeroPen_ * norm(Mc{1}(MZeros));
end  

if settings.mode_ == AltImplMode.EncourageDelay
    for t=1:Tc
        % formulating these constraints are slow, so output progress
        sprintf('Adding objectives for %d of %d matrices', t, Tc);
        
        BM = sys.B2 * Mc{t};        
        for i=1:sys.Nx
            for j=1:sys.Nx % note: this distance metric only for chain
                objective = objective + settings.fastCommPen_ * exp(abs(i-j)-t)*norm(Rc{t}(i,j));
                objective = objective + settings.fastCommPen_ * exp(abs(i-j)-t)*norm(BM(i,j));
            end
        end
    end
end

minimize(objective);
cvx_end

[Rc, Mc] = block_to_cell(Rcs, Mcs, Tc, sys);

if settings.mode_ == AltImplMode.StrictDelay
    for t=1:Tc % Even though these are constrained to be 0, set to 0
               % to avoid small numerical fluctuations
        Rc{t}(RZeros{t}) = 0;
        Mc{t}(MZeros{t}) = 0;
    end
end

if settings.mode_ == AltImplMode.StrictLocal
    for t=1:Tc % Even though these are constrained to be 0, set to 0
               % to avoid small numerical fluctuations
        Rc{t}(RZeros) = 0;
        Mc{t}(MZeros) = 0;
    end
end

% outputs
slsOuts_alt.R_           = Rc;
slsOuts_alt.M_           = Mc;
slsOuts_alt.clnorm_      = objective;
slsOuts_alt.solveStatus_ = cvx_status;
end


function slsOuts_alt = fai_analytic(sys, Tc, F1, F2)
% Note that this is meant to be used as a least-squares solver for 
% overconstrained systems; an exact non-optimized solution will be provided
% for a non-overconstrained system

% F2 * [Rcs; Mcs] == -F1;
RMc = F2 \ (-F1); % if overconstrained, least squares solution

% convert to cell form for easier calculation
Rcs      = RMc(1:sys.Nx * (Tc-1), :);
Mcs      = RMc(sys.Nx * (Tc-1) + 1:end, :);
[Rc, Mc] = block_to_cell(Rcs, Mcs, Tc, sys);

L1norm = norm([Rcs; Mcs], 1);
[slsOuts_alt.R_, slsOuts_alt.M_] = block_to_cell(Rcs, Mcs, Tc, sys);

slsOuts_alt.clnorm_      = L1norm;
slsOuts_alt.solveStatus_ = 'Analytic';
end