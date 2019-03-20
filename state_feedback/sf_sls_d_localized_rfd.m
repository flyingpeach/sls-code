function [R, M, acts] = sf_sls_d_localized_rfd(A,B,C,D,T,d,comms,ta,lambda,obj);

% System level synthesis with d-localizable constraints + actuator RFD
% Built around cvx
% Returns system response (R,M) as defined in Thm 1 of https://arxiv.org/pdf/1610.04815.pdf
Nx = size(A,1);
Nu = size(B,2);

[Rsupport,Msupport,count] = make_d_localized_constraints(A,B,T,d,comms,ta);

cvx_begin

variable X(count)
expression Rs(Nx,Nx,T)
expression Ms(Nu,Nx,T)

%populate decision variables
%locality constraints automatically enforced by limiting support of R and M
spot = 0;
for t = 1:T
    R{t} = Rs(:,:,t);
    supp = find(Rsupport{t});
    num = sum(sum(Rsupport{t}));
    R{t}(supp) = X(spot+1:spot+num);
    spot = spot + num;
    
    M{t} = Ms(:,:,t);
    supp = find(Msupport{t});
    num = sum(sum(Msupport{t}));
    M{t}(supp) = X(spot+1:spot+num);
    spot = spot + num;
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
    Mi = [];
    for t = 1:T
        Mi = [Mi,M{t}(i,:)];
    end
    act_penalty = act_penalty + norm(Mi,2);
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
    if  norm(vec(M{1}(i,:)),2)>1e-4
        acts = [acts;i];
    end
end
