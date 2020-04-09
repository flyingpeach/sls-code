% Algorithm I, System A

%% Setup
setup_system_a;
setup_sls_constr;
setup_loc_constr;

% Coupling weights and constraints
Q = eye(Nx);
S = diag(ones(Nu,1));

%% Controller synthesis

% ADMM variables
Phi = zeros(Nx*tFIR + Nu*(tFIR-1),Nx);
Psi = zeros(Nx*tFIR + Nu*(tFIR-1),Nx);
Lambda = zeros(Nx*tFIR + Nu*(tFIR-1),Nx);

x_t = x0;
for t = 1:tSim
    fprintf('Calculating time %d of %d\n', t, tSim); % display progress
    
    % ADMM parameters
    rho   = 5; % admm parameter
    eps_d = 1e-3; % convergence criterion for ||Psi(k+1) - Psi(k)||
    eps_p = 1e-4; % convergence criterion for ||Phi(k+1) - Psi(k+1)||
    maxIters = 5000;
    
    for iter=1:maxIters
        Psi_prev = Psi;
        
        %% Row-wise separability
        k = 0; % Separate the given matrices
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
        
        % Solve for each row            
        for i = 1:Nx
            clear ADMM_matrix
            ADMM_matrix = inv(2*x_t(s_r{i}(tFIR,:))*x_t(s_r{i}(tFIR,:))'+rho*eye(size(s_r{i},2)));
            Phi_loc{i} = rho*(Psi_loc_row{i}-Lambda_loc_row{i})*ADMM_matrix;
        end
        
        % Build the big matrix
        for i = 1:Nx
            Phi(r{i},s_r{i}(tFIR,:)) = Phi_loc{i};
        end
               
        %% Column-wise separability
        % Separate the given matrices
        for i = 1:Nx
            Phi_loc_col{i} = Phi(s_c{i},c{i});
            Lambda_loc_col{i} = Lambda(s_c{i},c{i});
        end
        
        % Solve for each column
        for i = 1:Nx
            clear AUX_matrix
            IZA_ZB_loc = IZA_ZB(:,s_c{i}); row_all_zeros = find(all(IZA_ZB_loc == 0,2)); keep_indices = setdiff(linspace(1,Nx*tFIR,Nx*tFIR),row_all_zeros);
            IZA_ZB_loc = IZA_ZB(keep_indices,s_c{i}); E1_loc = E1(keep_indices,c{i}); 
            AUX_matrix = IZA_ZB_loc'*pinv(IZA_ZB_loc*IZA_ZB_loc');
            Psi_loc{i} = (Phi_loc_col{i}+Lambda_loc_col{i})+AUX_matrix*(E1_loc-IZA_ZB_loc*(Phi_loc_col{i}+Lambda_loc_col{i}));
        end
         
        % Build the big matrix
        for i = 1:Nx
            Psi(s_c{i},c{i}) = Psi_loc{i};
        end
                     
        %% Lagrange multiplier
        Lambda = Lambda + Phi - Psi;
        
        % Check convergence locally 
        conv_criterion_failed = false;
        for sys = 1:Nx
            local_phi = Phi(r{sys},s_r{sys}(tFIR,:));
            local_psi = Psi(r{sys},s_r{sys}(tFIR,:));
            local_psi_prev = Psi_prev(r{sys},s_r{sys}(tFIR,:));

            local_diff_d = norm(local_psi-local_psi_prev,'fro');
            local_diff_p = norm(local_phi-local_psi,'fro');
            
            if local_diff_p > eps_p || local_diff_d > eps_d
                conv_criterion_failed = true;
                break; % if one fails, can stop checking the rest
            end
        end
        
        if ~conv_criterion_failed
            break; % convergence criterion passed, exit admm iterations
        end
    end

    if conv_criterion_failed
        fprintf('ADMM reached %d iters and did not converge\n', maxIters);
    end
    
    %% Dynamics
    % Compute the control action (in a localized way)
    u(:,t) = Phi(1+Nx*tFIR:Nx*tFIR+Nu,:)*x_t;
    
    % Simulate what the dynamics are given that action
    x(:,t+1) = Phi(1+Nx:2*Nx,:)*x_t; % Since there is no noise x_ref = x
    
    % Update the initial condition
    x_t = x(:,t+1);
    
end

%% Validation

x_VAL(:,1) = x0;
x_t = x0;

for k = 1:tSim
    
    clear LocalityR LocalityM
    
    Comms_Adj = abs(A)>0;
    LocalityR = Comms_Adj^(d-1)>0;
    
    count = 0;
    for t = 1:tFIR
        % Rsupport{t} = min(Comms_Adj^(floor(max(0,comms*(t-ta)))),LocalityR)>0;
        Rsupport{t} = LocalityR>0;
        Msupport{t} = (abs(B)'*Rsupport{t})>0;
        count = count + sum(sum(Rsupport{t}))+sum(sum(Msupport{t}));
    end
    
    cvx_begin
    cvx_precision low
    
    variable X(count)
    expression Rs(Nx,Nx,tFIR)
    expression Ms(Nu,Nx,tFIR)
    
    % Populate decision variables
    % Locality constraints automatically enforced by limiting support of R and M
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
        vect = vec([Q zeros(Nx,Nu); zeros(Nu,Nx) S]*[R{t};M{t}]*x_t);
        objective = objective + vect'*vect;
    end
    
    % Perform minimization
    minimize(objective)
    subject to
    % Achievability constraints
    R{1} == eye(Nx);
    for t= 1:tFIR-1
        R{t+1} == A*R{t} + B*M{t};
    end
    cvx_end
    
    %% Dynamics
    
    % Compute the control action
    u_VAL(:,k) = M{1}*x_t;
    
    % Simulate what the dynamics are given that action
    x_VAL(:,k+1) = R{2}*x_t; % Since there is no noise x_ref = x
    
    % Update the initial condition
    x_t = x_VAL(:,k+1); 
    
end

%% Cost 

obj=0;
for t =1:tSim
obj = obj + x(:,t)'*Q*x(:,t)+u(:,t)'*S*u(:,t);
end
obj = obj + x(:,t+1)'*Q*x(:,t+1);


obj_VAL=0;
for t =1:tSim
obj_VAL = obj_VAL + x_VAL(:,t)'*Q*x_VAL(:,t)+u_VAL(:,t)'*S*u_VAL(:,t);
end
obj_VAL = obj_VAL + x_VAL(:,t+1)'*Q*x_VAL(:,t+1);

obj-obj_VAL

% Output
header1 = 'Distributed MPC';
header2 = 'Centralized MPC!';
fprintf([ header1 ' ' header2 'r\n']);
fprintf('%f %f r\n', [obj obj_VAL]');

%% Plot

figure(1)
plot(1:tSim+1,x_VAL(1,:),'b',1:tSim+1,x(1,:),'*b',1:tSim+1,x_VAL(3,:),'g',1:tSim+1,x(3,:),'*g')
xlabel('$$Time$$','interpreter','latex','Fontsize', 16)
ylabel('$$\theta_{1},\ \theta_{2}$$','Interpreter','Latex','Fontsize', 16)
leg1 = legend('$$\theta_{1}\ Centralized\ MPC$$', '$$\theta_{1}\ Localized\ MPC\ using\ ADMM$$','$$\theta_{2}\ Centralized\ MPC$$', '$$\theta_{2}\ Localized\ MPC\ using\ ADMM$$');
set(leg1,'Interpreter','latex'); set(leg1, 'Fontsize', 10)
title('Subsystems 1 and 2')