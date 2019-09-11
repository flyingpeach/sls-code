function slsOuts_alt = find_alt_impl(sys, slsParams, slsOuts, Tc, approx)
% Find alternate implementation, returned in slsOuts_alt
% TODO: technically slsOuts_alt is a misnomer here as we didn't
%       get it from SLS but from postprocessing SLS
% Outputs
%    slsOuts_alt : contains clnorm and new R, M
% Inputs
%    sys         : LTISystem containing system matrices
%    slsParams   : SLSParams containing parameters
%    slsOuts     : contains info from SLS (original R, M)
%    approx      : put any value here if you want approximate solution
%                  (i.e. drop trailing zero constraints)
if nargin == 4
    approx    = false;
    statusTxt = 'Finding *exact* alt implementation';
else
    approx    = true;
    statusTxt = 'Finding *approx* alt implementation';
end

statusTxt = [char(10), statusTxt, sprintf(' with Tc=%d', Tc)];
disp(statusTxt);

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
    objective = objective + norm([sys.C1, sys.D12]*[Rc{t}; Mc{t}], 1);
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
if not(approx) % do not included these for approxmation
    for t=Tc+1:Tc+T
        for k=1:min(T, t)
             Rsum{t} = Rsum{t} + slsOuts.R_{k} * Dellc{t-k+1};
             Msum{t} = Msum{t} + slsOuts.M_{k} * Dellc{t-k+1};
        end
        Rsum{t} == 0;
        Msum{t} == 0;
    end
end

minimize(objective);
cvx_end

% outputs
slsOuts_alt.R_           = Rc;
slsOuts_alt.M_           = Mc;
slsOuts_alt.clnorm_      = objective;
slsOuts_alt.solveStatus_ = cvx_status;
