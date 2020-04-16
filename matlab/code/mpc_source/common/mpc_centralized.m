function [x, u, avgTime] = mpc_centralized(sys, d, Q, S, tFIR, tSim, x0, ...
                               varargin) %up, low/Ksmall, coupling
up=[]; low=[]; coupling=[];
try 
    up  = varargin{1};
end
try
    low = varargin{2};
end
try 
    coupling = varargin{3};
end
            
Nx = sys.Nx; Nu = sys.Nu; A = sys.A; B = sys.B2;

x      = zeros(Nx, tSim);
u      = zeros(Nu, tSim);
x(:,1) = x0;

totalTime = 0;

for t = 1:tSim
    fprintf('Validating time %d of %d\n', t, tSim); % display progress
    x_t = x(:,t);

    clear LocalityR LocalityM
    Comms_Adj = abs(A)>0;
    LocalityR = Comms_Adj^(d-1)>0;
    
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
    for k= 1:tFIR-1
        R{k+1} == A*R{k} + B*M{k};
    end
    
    if ~isempty(up) && ~isempty(low) 
        if isempty(coupling)
            % Bounding constraints specified
            for k = 1:tFIR
                R{k}*x_t <= up*ones(Nx,1);
                R{k}*x_t >= low*ones(Nx,1);
            end
        else
            % Coupling constraints specified
            Ksmall = low; % ULTRA HACKY
            for k = 3:tFIR
                Ksmall*R{k}*x_t <= up*ones(2*Nx,1);
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

avgTime = totalTime / (tSim - 1);

end