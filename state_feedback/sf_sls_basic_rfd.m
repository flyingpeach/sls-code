function [R, M, acts] = sf_sls_basic_rfd(A,B,C,D,T,lambda,obj);

% Basic version of system level synthesis state-feedback problem + actuator
% RFD
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

act_penalty = 0;
for i = 1:Nu
    act_penalty = act_penalty + norm(vec(Ms(i,:,:)),2);
end

minimize(objective + lambda*act_penalty)
subject to
    %achievability constraints
    R{1} == eye(Nx);

    for t= 1:T-1
        R{t+1} == A*R{t} + B*M{t};
    end
    R{T} == zeros(Nx,Nx);
cvx_end

acts = [];
for i=1:Nu
    if  norm(vec(Ms(i,:,1)),2)>1e-4
        acts = [acts;i];
    end
end

