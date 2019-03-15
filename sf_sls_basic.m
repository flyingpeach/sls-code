function [R, M, clnorm] = sf_sls_basic(A,B,C,D,T,obj);

% Basic version of system level synthesis state-feedback problem
% No locality constraints, only FIR
% Built around cvx
% Returns system response (R,M) as defined in Thm 1 of https://arxiv.org/pdf/1610.04815.pdf
Nx = size(A,1);
Nu = size(B,2);

cvx_begin
variable Rs(Nx,Nx,T)
variable Ms(Nu,Nx,T)

% populate decision variables
% this step isn't necessary but makes things easier to work with,
% especially when there are locality constraints, so including it here for
% compatibility reasons

for t = 1:T
    R{t} = Rs(:,:,t);
    M{t} = Ms(:,:,t);
end

objective = 0;

%set up objective function
switch obj
    case 'H2'
        objective = compute_H2(R,M,C,D,T);
    case 'Hinf'
        % todo
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
    %achievability constraints
    R{1} == eye(Nx);

    for t= 1:T-1
        R{t+1} == A*R{t} + B*M{t};
    end
    R{T} == zeros(Nx,Nx);
cvx_end

clnorm = objective;

