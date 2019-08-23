function [R, M, clnorm] = sf_sls_basic(A, B, C, D, TFIR, obj, xDesired)
% Basic system level synthesis with FIR only
% Returns 
%    R, M       : system response as defined in Thm 1 of https://arxiv.org/pdf/1610.04815.pdf
%    clnorm     : final (optimal value) of objective
% Inputs
%    A, B, C, D : system matrices with respect to regulated output 
%                 see (4.20) of https://arxiv.org/pdf/1904.01634.pdf
%    TFIR       : finite impulse response time
%    obj        : type of objective (i.e. 'H2', 'Hinf')
%    xDesired   : only used if obj is 'traj_track'; desired trajectory

Nx = size(A,1);
Nu = size(B,2);

cvx_begin
variable Rs(Nx,Nx,TFIR)
variable Ms(Nu,Nx,TFIR)

% populate decision variables

% this step isn't necessary but makes things easier to work with,
% especially when there are locality constraints, so including it here for
% compatibility reasons

for t = 1:TFIR
    R{t} = Rs(:,:,t);
    M{t} = Ms(:,:,t);
end

objective = 0;

%set up objective function
switch obj
    case 'traj_track'
        objective = traj_track(R, M, TFIR, Nx, Nu, xDesired);
    case 'H2'
        objective = compute_H2(R, M, C, D, TFIR);
    case 'Hinf'
        objective = compute_Hinf(R, M, C, D, TFIR);
    case 'L1'
        % todo
    case 'L1T'
        % todo
    otherwise
        % todo: throw a warning that only fidning a feasible solution
end
cvx_precision low
minimize(objective)
subject to
    % achievability constraints
    R{1} == eye(Nx);

    for t= 1:TFIR-1
        R{t+1} == A*R{t} + B*M{t};
    end
    R{TFIR} == zeros(Nx,Nx);
cvx_end

clnorm = objective;

