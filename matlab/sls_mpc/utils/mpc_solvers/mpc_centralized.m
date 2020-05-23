function [x, u, avgTime] = mpc_centralized(sys, x0, params)

fprintf('Centralized MPC\n')
params.sanity_check_cent();

Nx = sys.Nx; Nu = sys.Nu; A = sys.A; B = sys.B2;

locality = params.locality_;
tFIR     = params.tFIR_;   
tHorizon = params.tHorizon_;          

QSqrt = params.QSqrt_;
RSqrt = params.RSqrt_; % note this is not related to SLS R

x      = zeros(Nx, tHorizon);
u      = zeros(Nu, tHorizon);
x(:,1) = x0;

totalTime = 0;

% sparsity / support is locality-only, doesn't change with time
Comms_Adj = abs(A)>0;
RSupp     = Comms_Adj^(locality-1) > 0;
MSupp     = (abs(sys.B2)' * RSupp) > 0;
suppR     = find(RSupp);
suppM     = find(MSupp);
count     = tFIR * (length(suppR) + length(suppM));

for t = 1:tHorizon
    fprintf('Calculating time %d of %d\n', t, tHorizon); % display progress
    x_t = x(:,t);    
  
    cvx_begin quiet
    cvx_precision low
    
    expression Rs(Nx,Nx,tFIR)
    expression Ms(Nu,Nx,tFIR)
    variable RMSupp(count)
    
    spot = 0;
    for k = 1:tFIR
        R{k} = Rs(:,:,k); M{k} = Ms(:,:,k);

        R{k}(suppR) = RMSupp(spot+1:spot+length(suppR));
        spot        = spot + length(suppR);        
        M{k}(suppM) = RMSupp(spot+1:spot+length(suppM));
        spot        = spot + length(suppM);
    end
    
    % Set up objective function
    objective = 0;
    for k = 1:tFIR
        vect = vec([QSqrt zeros(Nx,Nu); zeros(Nu,Nx) RSqrt]*[R{k};M{k}]*x_t);
        objective = objective + vect'*vect;
    end

    tic;
    minimize(objective)
    subject to

    R{1} == eye(Nx); % Achievability constraints
    for k=1:tFIR-1
        R{k+1} == A*R{k} + B*M{k};
    end
    
    if params.has_state_cons()
        for k=2:tFIR % since we only care about later states
            if ~isinf(params.stateUB_)
                params.stateConsMtx_*R{k}*x_t <= params.stateUB_*ones(Nx, 1);
            end
            
            if ~isinf(params.stateLB_)
                params.stateConsMtx_*R{k}*x_t >= params.stateLB_*ones(Nx, 1);
            end
        end
    end

    cvx_end
    
    % Compute control + state
    u(:,t) = M{1}*x_t;
    x(:,t+1) = R{2}*x_t; % Since there is no noise x_ref = x 
    
    if k > 1
        totalTime = totalTime + toc;
    end
end

avgTime = totalTime / (tHorizon - 1);

end