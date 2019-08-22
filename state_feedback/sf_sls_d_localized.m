function [R, M, clnorm] = sf_sls_d_localized(A, B, C, D, TFIR, d, comms, ta, obj, xDesired)
% System level synthesis with d-localizable constraints
% Returns 
%    R, M       : system response as defined in Thm 1 of https://arxiv.org/pdf/1610.04815.pdf
%    clnorm     : final (optimal value) of objective

%    A, B, C, D : system matrices with respect to regulated output 
%                 see (4.20) of https://arxiv.org/pdf/1904.01634.pdf
%    TFIR       : finite impulse response time
%    d          : limit disturbance to d-hops 
%    comms      : communication speed
%    ta         : actuation delay
%    obj        : type of objective (i.e. 'H2', 'Hinf')
%    xDesired   : only used if obj is 'traj_track'; desired trajectory

Nx = size(A,1);
Nu = size(B,2);

[Rsupport,Msupport,count] = make_d_localized_constraints(A,B,TFIR,d,comms,ta);

cvx_begin
cvx_precision low
variable X(count)
expression Rs(Nx,Nx,TFIR)
expression Ms(Nu,Nx,TFIR)

% populate decision variables
% locality constraints automatically enforced by limiting support of R and M
spot = 0;
for t = 1:TFIR
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

