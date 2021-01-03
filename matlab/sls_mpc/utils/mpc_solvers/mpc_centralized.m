function [x, u, time] = mpc_centralized(sys, x0, params)

params.sanity_check_cent();

Nx = sys.Nx; Nu = sys.Nu; A = sys.A; B = sys.B2;

locality = params.locality_;
T     = params.tFIR_;   

QSqrt = params.QSqrt_;
RSqrt = params.RSqrt_; % note this is not related to SLS R

time = 0;

% sparsity / support is locality-only, doesn't change with time
Comms_Adj = abs(A)>0;
RSupp     = Comms_Adj^(locality-1) > 0;
MSupp     = (abs(sys.B2)' * RSupp) > 0;
suppR     = find(RSupp); % indices
suppM     = find(MSupp);
count     = length(suppR)*T + length(suppM)*(T-1);

cvx_begin quiet
cvx_precision low
    
expression Rs(Nx,Nx,T)
expression Ms(Nu,Nx,T-1)
variable RMSupp(count)

R = cell(T, 1);
M = cell(T-1, 1);

spot = 0;
for k = 1:T
    R{k}        = Rs(:,:,k);
    R{k}(suppR) = RMSupp(spot+1:spot+length(suppR));
    spot        = spot + length(suppR);
    if k < T
        M{k}        = Ms(:,:,k);
        M{k}(suppM) = RMSupp(spot+1:spot+length(suppM));
        spot        = spot + length(suppM);
    end
end
    
% Set up objective function
objective = 0;
for k = 1:T
    if k < T
        M_ = M{k};
    else
        M_ = zeros(Nu,Nx);
    end
    vect = vec([QSqrt zeros(Nx,Nu); zeros(Nu,Nx) RSqrt]*[R{k};M_]*x0);
    objective = objective + vect'*vect;
end

tic;
minimize(objective)
subject to

R{1} == eye(Nx); % Achievability constraints
for k=1:T-1
    R{k+1} == A*R{k} + B*M{k};
end
    
if params.has_state_cons()
    for k=1:T
        if ~isinf(params.stateUB_)
            params.stateConsMtx_*R{k}*x0 <= params.stateUB_*ones(Nx, 1);
        end
            
        if ~isinf(params.stateLB_)
            params.stateConsMtx_*R{k}*x0 >= params.stateLB_*ones(Nx, 1);
        end
    end
end
    
if params.has_input_cons()
    for k=1:T-1
        if ~isinf(params.inputUB_)
            params.inputConsMtx_*M{k}*x0 <= params.inputUB_*ones(Nu, 1);
        end
            
        if ~isinf(params.inputLB_)
            params.inputConsMtx_*M{k}*x0 >= params.inputLB_*ones(Nu, 1);
        end
    end        
end

time = time + toc;
cvx_end
    
% Compute control + state
u = M{1}*x0;
x = R{2}*x0; % Since there is no noise x_ref = x 

end