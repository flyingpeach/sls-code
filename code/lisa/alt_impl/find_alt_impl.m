function slsOuts_alt = find_alt_impl(sys, slsParams, slsOuts, Tc)
% Find alternate implementation, returned in slsOuts_alt
% TODO: technically slsOuts_alt is a misnomer here as we didn't
%       get it from SLS but from postprocessing SLS
%    slsOuts_alt : contains clnorm and new R, M
% Inputs
%    sys         : LTISystem containing system matrices
%    slsParams   : SLSParams containing parameters
%    slsOuts     : contains info from SLS (original R, M)

statusTxt = sprintf('Finding alternate CL implementation with Tc=%d', Tc);
statusTxt = [char(10), statusTxt];
disp(statusTxt);

cvx_begin
cvx_precision low

TMax = max(Tc, slsParams.tFIR_);

variable Rcs(sys.Nx, sys.Nx, Tc)
variable Mcs(sys.Nu, sys.Nx, Tc)
expression Dellcs(sys.Nx, sys.Nx, Tc)
expression Rcsums(sys.Nx, sys.Nx, Tc)
expression Mcsums(sys.Nu, sys.Nx, Tc)

% populate decision variables
objective = 0;
for t = 1:Tc
    Rc{t}     = Rcs(:,:,t);
    Mc{t}     = Mcs(:,:,t);
    Dellc{t}  = Dellcs(:,:,t);
    Rcsum{t}  = Rcsums(:,:,t);
    Mcsum{t}  = Mcsums(:,:,t);
    
    % L1 norm to enforce sparsity
    objective = objective + norm([sys.C1, sys.D12]*[Rc{t}; Mc{t}], 1);
end

% Not strictly necessary but makes everything easier
% Equivalent to Rc{1} == eye(sys.Nx)
% Should not reduce dimensionality of feasible Rc/Mc space
Dellc{1} = eye(sys.Nx);

for t=2:Tc-1
    Dellc{t} = Rc{t+1} - sys.A*Rc{t} - sys.B2*Mc{t};
end
Dellc{Tc} = -sys.A*Rc{Tc} - sys.B2*Mc{Tc};
for t=Tc+1:TMax
    Dellc{t} = zeros(sys.Nx, sys.Nx);
end

for t=1:Tc
    for k=1:min([slsParams.tFIR_, t]) % convolve
         Rcsum{t} = Rcsum{t} + slsOuts.R_{k} * Dellc{t-k+1};
         Mcsum{t} = Mcsum{t} + slsOuts.M_{k} * Dellc{t-k+1};
    end
    Rc{t} == Rcsum{t};
    Mc{t} == Mcsum{t};
end

minimize(objective);
cvx_end

% outputs
slsOuts_alt.R_      = Rc;
slsOuts_alt.M_      = Mc;
slsOuts_alt.clnorm_ = objective;
