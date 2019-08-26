function [R, M, clnorm] = sf_sls_basic(sys, params)
% Basic system level synthesis with FIR only
% Returns 
%    R, M   : system response as defined in Thm 1 of https://arxiv.org/pdf/1610.04815.pdf
%    clnorm : final (optimal value) of objective
% Inputs
%    sys    : an LTISystem
%    params : SLSParams containing parameters

cvx_begin
variable Rs(sys.Nx, sys.Nx, params.tFIR_)
variable Ms(sys.Nu, sys.Nx, params.tFIR_)

% populate decision variables
% this step isn't necessary but makes things easier to work with
for t = 1:params.tFIR_
    R{t} = Rs(:,:,t);
    M{t} = Ms(:,:,t);
end

objective = get_objective(sys, params, R, M);

% solve optimization problem
cvx_precision low
minimize(objective)
subject to
    % achievability constraints
    R{1} == eye(sys.Nx);

    for t= 1:params.tFIR_-1
        R{t+1} == sys.A*R{t} + sys.B2*M{t};
    end
    R{params.tFIR_} == zeros(sys.Nx, sys.Nx);
cvx_end

clnorm = objective;

