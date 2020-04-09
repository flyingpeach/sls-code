function [x, u] = mpc_centralized(Nx, Nu, A, B, d, Q, S, tFIR, tSim, x0)

x      = zeros(Nx, tSim);
u      = zeros(Nu, tSim);
x(:,1) = x0;

for k = 1:tSim
    fprintf('Validating time %d of %d\n', k, tSim); % display progress
    x_k = x(:,k);

    clear LocalityR LocalityM
    Comms_Adj = abs(A)>0;
    LocalityR = Comms_Adj^(d-1)>0;
    
    count = 0;
    for t = 1:tFIR
        Rsupport{t} = LocalityR>0;
        Msupport{t} = (abs(B)'*Rsupport{t})>0;
        count = count + sum(sum(Rsupport{t}))+sum(sum(Msupport{t}));
    end
    
    cvx_begin quiet
    cvx_precision low
    
    variable X(count)
    expression Rs(Nx,Nx,tFIR)
    expression Ms(Nu,Nx,tFIR)
    
    R = cell(1, tFIR);
    M = cell(1, tFIR);
    
    spot = 0;
    for t = 1:tFIR
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
    
    % Set up objective function
    objective = 0;
    for t = 1:tFIR
        vect = vec([Q zeros(Nx,Nu); zeros(Nu,Nx) S]*[R{t};M{t}]*x_k);
        objective = objective + vect'*vect;
    end

    minimize(objective)
    subject to

    R{1} == eye(Nx); % Achievability constraints
    for t= 1:tFIR-1
        R{t+1} == A*R{t} + B*M{t};
    end
    cvx_end
    
    % Compute control + state
    u(:,k) = M{1}*x_k;
    x(:,k+1) = R{2}*x_k; % Since there is no noise x_ref = x 
end

end