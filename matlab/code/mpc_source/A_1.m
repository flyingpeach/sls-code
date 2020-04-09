% Algorithm I, System A

%% Setup
setup_system_a;
setup_sls_constr;
setup_loc_constr;

% Weights
Q = eye(Nx);
S = diag(ones(Nu,1));

%% Distributed MPC

% ADMM variables
Phi    = zeros(Nx*tFIR + Nu*(tFIR-1),Nx);
Psi    = zeros(Nx*tFIR + Nu*(tFIR-1),Nx);
Lambda = zeros(Nx*tFIR + Nu*(tFIR-1),Nx);

% Keep track of state + control
x       = zeros(Nx, tSim);
u       = zeros(Nu, tSim);
x(:, 1) = x0;

for t = 1:tSim
    fprintf('Calculating time %d of %d\n', t, tSim); % display progress
    x_t = x(:,t); % update initial condition
    
    % ADMM parameters
    rho   = 5; % admm parameter
    eps_d = 1e-3; % convergence criterion for ||Psi(k+1) - Psi(k)||
    eps_p = 1e-4; % convergence criterion for ||Phi(k+1) - Psi(k+1)||
    maxIters = 5000;
    
    for iter=1:maxIters
        Psi_prev = Psi;
        
        Psi_loc_row    = cell(1, Nx);
        Lambda_loc_row = cell(1, Nx);
        
        % Separate into rows
        k = 0;
        for i = 1:Nx
            if mod(i, Nx/Nu) == 0
                 k = k+1;
                for j = 1:tFIR+(tFIR-1)
                    if j<=tFIR
                        Psi_loc_row{i} = Psi(r{i},s_r{i}(j,1:max(length(find(LocalityR{j}(i,:))))));
                        Lambda_loc_row{i} = Lambda(r{i},s_r{i}(j,1:max(length(find(LocalityR{j}(i,:))))));
                    else
                        Psi_loc_row{i} = Psi(r{i},s_r{i}(j,1:max(length(find(LocalityM{j-tFIR}(k,:))))));
                        Lambda_loc_row{i} = Lambda(r{i},s_r{i}(j,1:max(length(find(LocalityM{j-tFIR}(k,:))))));
                    end
                end
            else
                for j = 1:tFIR
                        Psi_loc_row{i} = Psi(r{i},s_r{i}(j,1:max(length(find(LocalityR{j}(i,:))))));
                        Lambda_loc_row{i} = Lambda(r{i},s_r{i}(j,1:max(length(find(LocalityR{j}(i,:))))));
                end
            end
        end
        
        Phi_loc = cell(1, Nx);
        
        % Solve for each row            
        for i = 1:Nx
            clear ADMM_matrix
            ADMM_matrix = inv(2*x_t(s_r{i}(tFIR,:))*x_t(s_r{i}(tFIR,:))'+rho*eye(size(s_r{i},2)));
            Phi_loc{i} = rho*(Psi_loc_row{i}-Lambda_loc_row{i})*ADMM_matrix;
        end
        
        % Build Phi matrix
        for i = 1:Nx
            Phi(r{i},s_r{i}(tFIR,:)) = Phi_loc{i};
        end
        
        Phi_loc_col    = cell(1, Nx);
        Lambda_loc_col = cell(1, Nx);               
        % Separate into columns
        for i = 1:Nx
            Phi_loc_col{i} = Phi(s_c{i},c{i});
            Lambda_loc_col{i} = Lambda(s_c{i},c{i});
        end
        
        Psi_loc = cell(1, Nx);
        
        % Solve for each column
        for i = 1:Nx
            clear AUX_matrix
            IZA_ZB_loc = IZA_ZB(:,s_c{i}); row_all_zeros = find(all(IZA_ZB_loc == 0,2)); keep_indices = setdiff(linspace(1,Nx*tFIR,Nx*tFIR),row_all_zeros);
            IZA_ZB_loc = IZA_ZB(keep_indices,s_c{i}); E1_loc = E1(keep_indices,c{i}); 
            AUX_matrix = IZA_ZB_loc'*pinv(IZA_ZB_loc*IZA_ZB_loc');
            Psi_loc{i} = (Phi_loc_col{i}+Lambda_loc_col{i})+AUX_matrix*(E1_loc-IZA_ZB_loc*(Phi_loc_col{i}+Lambda_loc_col{i}));
        end
         
        % Build Phi matrix
        for i = 1:Nx
            Psi(s_c{i},c{i}) = Psi_loc{i};
        end
                     
        % Lagrange multiplier
        Lambda = Lambda + Phi - Psi;
        
        % Check convergence locally 
        check_conv;
        if ~criterion_failed
            break; % convergence criterion passed, exit admm iterations
        end
    end

    if criterion_failed
        fprintf('ADMM reached %d iters and did not converge\n', maxIters);
    end
    
    % Compute control + state
    u(:,t) = Phi(1+Nx*tFIR:Nx*tFIR+Nu,:)*x_t;
    x(:,t+1) = Phi(1+Nx:2*Nx,:)*x_t; % since no noise, x_ref = x
end

%% Centralized MPC (for validation + comparison)

xVal      = zeros(Nx, tSim);
uVal      = zeros(Nu, tSim);
xVal(:,1) = x0;

for k = 1:tSim
    fprintf('Validating time %d of %d\n', k, tSim); % display progress
    x_k = xVal(:,k);

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
    uVal(:,k) = M{1}*x_k;
    xVal(:,k+1) = R{2}*x_k; % Since there is no noise x_ref = x 
end

%% Calculate costs + plot 

obj=0;
for t=1:tSim
    obj = obj + x(:,t)'*Q*x(:,t)+u(:,t)'*S*u(:,t);
end
obj = obj + x(:,t+1)'*Q*x(:,t+1);

objVal=0;
for t=1:tSim
    objVal = objVal + xVal(:,t)'*Q*xVal(:,t)+uVal(:,t)'*S*uVal(:,t);
end
objVal = objVal + xVal(:,t+1)'*Q*xVal(:,t+1);

% Output costs
fprintf('Distributed cost: %f\n', obj);
fprintf('Centralized cost: %f\n', objVal);

figure(1)
plot(1:tSim+1,xVal(1,:),'b',1:tSim+1,x(1,:),'*b',1:tSim+1,xVal(3,:),'g',1:tSim+1,x(3,:),'*g')
xlabel('$$Time$$','interpreter','latex','Fontsize', 10)
ylabel('$$\theta_{1},\ \theta_{2}$$','Interpreter','Latex','Fontsize', 10)
leg1 = legend('$$\theta_{1}\ Centralized\ MPC$$', '$$\theta_{1}\ Localized\ MPC\ using\ ADMM$$','$$\theta_{2}\ Centralized\ MPC$$', '$$\theta_{2}\ Localized\ MPC\ using\ ADMM$$');
set(leg1,'Interpreter','latex'); set(leg1, 'Fontsize', 8)
title('Subsystems 1 and 2')