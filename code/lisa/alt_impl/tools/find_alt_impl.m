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

if settings.mode_ ~= AltImplMode.ImplicitOpt
    F  = get_F(sys, slsParams, slsOuts, Tc);
    F1 = F(:, 1:sys.Nx);
    F2 = F(:,sys.Nx+1:end);
end
    
switch settings.mode_
    case AltImplMode.ImplicitOpt
        slsOuts_alt = fai_implicit(sys, slsParams, slsOuts, Tc);
    case AltImplMode.ExplicitOpt
        leaky       = false;
        slsOuts_alt = fai_explicit(sys, Tc, F1, F2, leaky);
    case AltImplMode.Analytic
        slsOuts_alt = fai_analytic(sys, Tc, F1, F2);
    case AltImplMode.ApproxDrop
        numDrop     = settings.relaxPct_ * size(F, 1);
        [F1, F2]    = drop_constr(F1, F2, numDrop);
        leaky       = false;
        slsOuts_alt = fai_explicit(sys, Tc, F1, F2, leaky);
    case AltImplMode.ApproxLeaky
        leaky       = true;
        slsOuts_alt =  fai_explicit(sys, Tc, F1, F2, leaky, settings.clDiffPen_);
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


function [F1, F2] = drop_constr(F1, F2, numDrop)
% drop rows of F1 (and corresponding rows in F2) with lowest row norm
F1rownorms = sqrt(sum(F1.^2, 2));
[~, idx]   = sort(F1rownorms);

keepRows   = 1:size(F1, 1); % indices of rows to keep

for i=1:numDrop
    keepRows(keepRows==idx(i)) = []; % get rid of that row
end

F1 = F1(keepRows, :);
F2 = F2(keepRows, :);
end


function slsOuts_alt = fai_implicit(sys, slsParams, slsOuts, Tc)
T = slsParams.tFIR_;

cvx_begin
cvx_solver sdpt3
cvx_precision low

variable Rcs(sys.Nx, sys.Nx, Tc)
variable Mcs(sys.Nu, sys.Nx, Tc)
expression Dellcs(sys.Nx, sys.Nx, Tc+1)
expression Rsums(sys.Nx, sys.Nx, Tc+T)
expression Msums(sys.Nu, sys.Nx, Tc+T)

% populate decision variables
objective = 0;
for t = 1:Tc
    Rc{t}    = Rcs(:,:,t);
    Mc{t}    = Mcs(:,:,t);
    Dellc{t} = Dellcs(:,:,t);
    Rsum{t}  = Rsums(:,:,t);
    Msum{t}  = Msums(:,:,t);
    
    % L1 norm to enforce sparsity
    objective = objective + norm([Rc{t}; Mc{t}], 1);
end
Dellc{Tc+1} = Dellcs(:,:,Tc+1);
for t=Tc+1:Tc+T
    Rsum{t}  = Rsums(:,:,t);
    Msum{t}  = Msums(:,:,t);
end

% calculate Dellc
Dellc{1} = eye(sys.Nx); % will enforce Rc{1} == eye(sys.Nx)
for t=2:Tc
    Dellc{t} = Rc{t} - sys.A*Rc{t-1} - sys.B2*Mc{t-1};
end
Dellc{Tc+1} = -sys.A*Rc{Tc} - sys.B2*Mc{Tc};
for t=Tc+2:Tc+T
    Dellc{t} = zeros(sys.Nx, sys.Nx); % zero padding for ease of calculation
end

% enforce LHS = RHS constraints
for t=1:Tc
    for k=1:min(T, t) % convolve
         Rsum{t} = Rsum{t} + slsOuts.R_{k} * Dellc{t-k+1};
         Msum{t} = Msum{t} + slsOuts.M_{k} * Dellc{t-k+1};
    end
    Rc{t} == Rsum{t};
    Mc{t} == Msum{t};
end
for t=Tc+1:Tc+T % note: in initial prototypes, dropped these constr for approximation
    for k=1:min(T, t)
         Rsum{t} = Rsum{t} + slsOuts.R_{k} * Dellc{t-k+1};
         Msum{t} = Msum{t} + slsOuts.M_{k} * Dellc{t-k+1};
    end
        Rsum{t} == 0;
        Msum{t} == 0;
end

minimize(objective);
cvx_end

% outputs
slsOuts_alt.R_           = Rc;
slsOuts_alt.M_           = Mc;
slsOuts_alt.clnorm_      = objective;
slsOuts_alt.solveStatus_ = cvx_status;
end


function slsOuts_alt = fai_explicit(sys, Tc, F1, F2, leaky, clDiffPen)
cvx_begin
cvx_solver sdpt3
cvx_precision low

variable Rcs(sys.Nx * (Tc-1), sys.Nx) % Rc{1} == 1 so it's not a variable
variable Mcs(sys.Nu * Tc, sys.Nx)
expression Delta(Tc * (sys.Nx + sys.Nu) - sys.Nx, sys.Nx) 
    
objective = norm([Rcs; Mcs], 1);

if leaky
    Delta = F2 * [Rcs; Mcs] + F1;    
    clDiff = norm(Delta);
    objective = objective + clDiffPen * clDiff;
else % exact solution
    F2 * [Rcs; Mcs] == -F1;
end

minimize(objective);
cvx_end

% outputs
[slsOuts_alt.R_, slsOuts_alt.M_] = block_to_cell(Rcs, Mcs, Tc, sys);

% want original L1norm
if leaky
    objective = objective - clDiffPen * clDiff;
end
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