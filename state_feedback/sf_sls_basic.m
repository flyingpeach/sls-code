function [R, M, clnorm] = sf_sls_basic(sys, TFIR, obj, xDes)
% Basic system level synthesis with FIR only
% Returns 
%    R, M   : system response as defined in Thm 1 of https://arxiv.org/pdf/1610.04815.pdf
%    clnorm : final (optimal value) of objective
% Inputs
%    sys    : an LTISystem
%    TFIR   : finite impulse response time
%    obj    : type of objective (i.e. 'H2', 'Hinf')
%    xDes   : only used if obj is 'traj_track'; desired trajectory

cvx_begin
variable Rs(sys.Nx, sys.Nx, TFIR)
variable Ms(sys.Nu, sys.Nx, TFIR)

% populate decision variables
% this step isn't necessary but makes things easier to work with
for t = 1:TFIR
    R{t} = Rs(:,:,t);
    M{t} = Ms(:,:,t);
end

objective = 0;

% set up objective function
switch obj
    case 'traj_track'
        objective = traj_track(R, M, TFIR, sys.Nx, sys.Nu, xDes);
    case 'H2'
        objective = compute_H2(R, M, sys.C1, sys.D12, TFIR);
    case 'Hinf'
        objective = compute_Hinf(R, M, sys.C1, sys.D12, TFIR);
    case 'L1'
        % todo
    case 'L1T'
        % todo
    otherwise
        % todo: throw a warning that only fidning a feasible solution
end

% solve optimization problem
cvx_precision low
minimize(objective)
subject to
    % achievability constraints
    R{1} == eye(sys.Nx);

    for t= 1:TFIR-1
        R{t+1} == sys.A*R{t} + sys.B2*M{t};
    end
    R{TFIR} == zeros(sys.Nx, sys.Nx);
cvx_end

clnorm = objective;

