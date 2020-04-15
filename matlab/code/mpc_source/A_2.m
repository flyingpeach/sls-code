% Algorithm II, System A

%% Setup
setup_system_a;

% Weights
Q = zeros(Nx,Nx);
for i = 1:2:Nx
    Q(i,i) = 1;
    if i > 1
        Q(i,i-2) = -1/2;
    end
    if i < Nx-1
        Q(i,i+2) = -1/2;
    end
    if i < Nx
        Q(i+1,i+1) = .01;
    end
end
S = eye(Nu);

% Build cost matrix
C = [];
for t = 0:tFIR-1
    C = blkdiag(C, Q);
end
for t = 0:tFIR-2
    C = blkdiag(C, S);
end    

% Coupling
indices = cell(1, length(C));
for i = 1:length(C)
    for j = 1:length(C)
        if C(i,j) ~= 0
            indices{i} = [indices{i} j];
        end
    end
end
%% ADMM

%% Syntheize the controller

x(:,1) = x0;
xi = x0;

% Warm-start
Phi = zeros(Nx*tFIR + Nu*(tFIR-1),Nx);
Psi = zeros(Nx*tFIR + Nu*(tFIR-1),Nx);
Lambda = zeros(Nx*tFIR + Nu*(tFIR-1),Nx);
for i = 1:Nx*tFIR+Nu*(tFIR-1)
    if isempty(indices{i}) == 0
        Z_admm{i} = 0;
        for j = indices{i}
            Y{i}{j} = 0;
        end
    else
    end
end

for t = 1:tSim
    t
    
Psi_prev = ones(Nx*tFIR + Nu*(tFIR-1),Nx); % Just so the while doesn't break

rho = 1;

count = 0;conv = [1];
while norm(conv) ~= 0 
    
Psi_prev = Psi;

        %% Row-wise separability
        
        % Separate the given matrices
k = 0;
for sys = 1:Nx
    if mod(sys, Nx/Nu) == 0
        k = k+1;
    end
    for i = r{sys}
        j = find(r{sys}==i);
        if j<=tFIR
            Psi_loc_row{i} = Psi(i,s_r{sys}(j,1:max(length(find(LocalityR{j}(sys,:))))));
            Lambda_loc_row{i} = Lambda(i,s_r{sys}(j,1:max(length(find(LocalityR{j}(sys,:))))));
        else
            Psi_loc_row{i} = Psi(i,s_r{sys}(j,1:max(length(find(LocalityM{j-tFIR}(k,:))))));
            Lambda_loc_row{i} = Lambda(i,s_r{sys}(j,1:max(length(find(LocalityM{j-tFIR}(k,:))))));
        end
    end
end


%% ADMM
counter = 0; conv1 = 1; conv2 = 1;

mu = 1;

while norm(conv1)>10^(-3) || norm(conv2)>10^(-4)
    
    
    %% Phi & X
    
    for sys = 1:Nx
        
        for i = r{sys}
            
            n = max(length(s_r{sys}(find(r{sys}==i),:)));
            
            % Define parameters
            m_j = max(length(indices{i}));
            C_proxi = C(i,indices{i});
            M1 = [eye(n) zeros(n,m_j-1)];
            Id = eye(m_j-1);
            xi_i = xi(s_r{sys}(find(r{sys}==i),:));
            
            i_new = find(indices{i} == i);
            
            M2 = [zeros(i_new-1,n) Id(1:i_new-1,:); xi_i' zeros(1,m_j-1); zeros(m_j-i_new,n) Id(i_new:end,:)];
            a = Psi_loc_row{i}-Lambda_loc_row{i};
            
            Mj_sum = 0; Mjb_sum = 0;
            for j = 1:m_j
                
                if j < i_new
                    M{j} = zeros(1,m_j+n-1); M{j}(n+j) = 1;
                elseif j == i_new
                    M{j} = [xi_i' zeros(1,m_j-1)];
                elseif j>i_new
                    M{j} = zeros(1,m_j+n-1); M{j}(n+j-1) = 1;
                end
                
                Mj_sum = Mj_sum + (M{j}'*M{j});
                
                k = indices{i}(j);
                b = Z_admm{k}-Y{i}{k};
                Mjb_sum = Mjb_sum + M{j}'*b;
            end
            
            W = pinv(2*((C_proxi*M2)'*C_proxi*M2)+rho*(M1'*M1)+mu*Mj_sum)*(rho*M1'*a'+mu*Mjb_sum);
            
            
            Phi_loc{i} = (M1*W)';
            
            X_i = M2*W;
            X{i} = zeros(Nx*tFIR+Nu*(tFIR-1),1); X{i}(indices{i}) = X_i;
            
        end
    end
    
    %% Z
    
    Z_old = Z_admm;
     
    for sys = 1:Nx
        for j = r{sys}
            Z_admm{j} = 0;
            for i = indices{j}
                Z_admm{j} = Z_admm{j} + (X{i}(j)+Y{i}{j})/length(indices{j});
            end
        end
    end
    
    %% Lagrange multiplier
    
    for sys = 1:Nx
        for i = r{sys}
            for j = indices{i}
                Y{i}{j} = Y{i}{j} + X{i}(j) - Z_admm{j};
            end
        end
    end
    
    %% Convergence criterium
    
    % Primal residue
    for sys = 1:Nx
        for i = r{sys}
            average = zeros(Nx*tFIR+Nu*(tFIR-1),1);
            for j = indices{i}
                average(j) = Z_admm{j};
            end
            primal (i) = norm(X{i}-average);
        end
    end
    conv1 = max(primal); 

    % Dual residue 
    for sys = 1:Nx
        for j = r{sys}
            dual(j) = norm(Z_admm{j}-Z_old{j});
        end
    end
    conv2 = max(dual); 
    
    counter = counter +1;
    if counter > 200
        break
        disp('Consensus ADMM did not converge')
    end

end

% Build the big matrix
for sys = 1:Nx
    for i = r{sys}
        Phi(i,s_r{sys}(find(r{sys}==i),:)) = Phi_loc{i};
    end
end

    %% Column-wise separability
        % Separate the given matrices
        for i = 1:Nx
            Phi_loc_col{i} = Phi(s_c{i},c{i});
            Lambda_loc_col{i} = Lambda(s_c{i},c{i});
        end

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

    %% Convergence
        % Local convergence criterium
        conv = [0];
        
        for sys = 1:Nx
            local_phi = Phi(r{sys},s_r{sys}(tFIR,:));
            local_psi = Psi(r{sys},s_r{sys}(tFIR,:));
            local_psi_prev = Psi_prev(r{sys},s_r{sys}(tFIR,:));

            local_conv1 = norm(local_phi-local_psi,'fro');
            local_conv2 = norm(local_psi-local_psi_prev,'fro');
            
            if local_conv1 > 10^(-4) || local_conv2 > 10^(-3)
                 conv = [conv 1];
            end
        end
        
        % Number of iterations until convergence
        count = count + 1;
        if count > 5000
            disp ('ADMM did not converge')
            break
        end

end

%% Dynamics
    
    % Compute the control action (in a localized way)
    u(:,t) = Phi(1+Nx*tFIR:Nx*tFIR+Nu,:)*xi;
    
    % Simulate what the dynamics are given that action
    x(:,t+1) = Phi(1+Nx:2*Nx,:)*xi; % Since there is no noise x_ref = x
    
    % Update the initial condition
    xi = x(:,t+1);
    
end

%% Centralized MPC (for validation + comparison)
[xVal, uVal, ~] = mpc_centralized(Nx, Nu, A, B, d, Q, S, tFIR, tSim, x0);

%% Calculate costs + plot 
obj    = get_cost_fn(Q, S, tSim, x, u);
objVal = get_cost_fn(Q, S, tSim, xVal, uVal);

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