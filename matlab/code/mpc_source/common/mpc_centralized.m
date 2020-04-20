function [x, u, avgTime] = mpc_centralized(sys, x0, params, Q, S)

Nx = sys.Nx; Nu = sys.Nu; A = sys.A; B = sys.B2;

locality = params.locality_;
tFIR     = params.tFIR_;   
tHorizon = params.tHorizon_;          

x      = zeros(Nx, tHorizon);
u      = zeros(Nu, tHorizon);
x(:,1) = x0;

totalTime = 0;

for t = 1:tHorizon
    fprintf('Validating time %d of %d\n', t, tHorizon); % display progress
    x_t = x(:,t);

    clear LocalityR LocalityM
    Comms_Adj = abs(A)>0;
    LocalityR = Comms_Adj^(locality-1)>0;
    
    count = 0;
    for k = 1:tFIR
        Rsupport{k} = LocalityR>0;
        Msupport{k} = (abs(B)'*Rsupport{k})>0;
        count = count + sum(sum(Rsupport{k}))+sum(sum(Msupport{k}));
    end
    
    cvx_begin quiet
    cvx_precision low
    
    variable X(count)
    expression Rs(Nx,Nx,tFIR)
    expression Ms(Nu,Nx,tFIR)
    
    R = cell(1, tFIR);
    M = cell(1, tFIR);
    
    spot = 0;
    for k = 1:tFIR
        R{k} = Rs(:,:,k);
        supp = find(Rsupport{k});
        num = sum(sum(Rsupport{k}));
        R{k}(supp) = X(spot+1:spot+num);
        spot = spot + num;
        
        M{k} = Ms(:,:,k);
        supp = find(Msupport{k});
        num = sum(sum(Msupport{k}));
        M{k}(supp) = X(spot+1:spot+num);
        spot = spot + num;
    end
    
    % Set up objective function
    objective = 0;
    for k = 1:tFIR
        vect = vec([Q zeros(Nx,Nu); zeros(Nu,Nx) S]*[R{k};M{k}]*x_t);
        objective = objective + vect'*vect;
    end

    tic;
    minimize(objective)
    subject to

    R{1} == eye(Nx); % Achievability constraints
    for k=1:tFIR-1
        R{k+1} == A*R{k} + B*M{k};
    end
    
    if ~isempty(params.state_upperbnd_)
        for k=1:tFIR
            R{k}*x_t <= params.state_upperbnd_*ones(Nx,1);
        end
    end    
    if ~isempty(params.state_lowerbnd)
        for k=1:tFIR
            R{k}*x_t >= params.state_lowerbnd_*ones(Nx,1);
        end
    end

    if ~isempty(params.couplingMtx_)
        for k = 3:tFIR % TODO: why 3?
            params.couplingMtx_*R{k}*x_t <= params.state_upperbnd*ones(2*Nx,1);
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