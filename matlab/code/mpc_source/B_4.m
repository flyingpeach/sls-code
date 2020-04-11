% Algorithm II, System B
clear all; close all; clc;

localities = [3 5 7 10];
numLocs    = length(localities);
times      = zeros(1, numLocs);
timeCents  = zeros(1, numLocs);
iters      = zeros(1, numLocs);

numPendula = 10;
for index=1:numLocs
    locality = localities(index);
    clear x u x_VAL u_VAL LocalityR LocalityM
    
    setup_system_b;
    setup_sls_constr;
    setup_loc_constr;

%% Coupling weights and constraints

% Weights
Q = diag(ones(Nx,1))+diag(-1/2*ones(Nx-2,1),2)+diag(-1/2*ones(Nx-2,1),-2);
S = diag(ones(Nu,1));

% Build the big cost matrix
C = [];
for t = 0:tFIR-1
    C =  [C; zeros(Nx,t*Nx) Q zeros(Nx,(tFIR-t-1)*Nx+(tFIR-1)*Nu)];
end
for t = 0:tFIR-2
    C =  [C;zeros(Nu,tFIR*Nx+t*Nu) S zeros(Nu,(tFIR-t-2)*Nu)];
end

% Constraints
for i = 1:2:2*(Nx-1)
    K1(i,i) = 1; K1(i,i+2) = -1;
    K1(i+1,i) = -1; K1(i+1,i+2) = 1;
end
K1 = K1(1:Nx,1:Nx); K1(Nx-1:Nx,:) = zeros(2,Nx);
  
Ksmall = zeros(2*Nx,Nx); j = 0;
for i = 1:2*Nx
    if mod(i,4) == 1 | mod(i,4) == 2
        j = j + 1;
        Ksmall(i,:) = K1(j,:);
    else
    end
end

% Build the big constraint matrix
K = [zeros(size(Ksmall)) zeros(2*Nx,(tFIR-1)*Nx)];
K =  [K; zeros(2*Nx,Nx) zeros(size(Ksmall)) zeros(2*Nx,(tFIR-2)*Nx)];
for t = 2:tFIR-1
    K =  [K; zeros(2*Nx,t*Nx) Ksmall zeros(2*Nx,(tFIR-t-1)*Nx)];
end

% Upper bound
up = .5;

% Coupling (since the coupling from cost is more extensive than the coupling from the constraints this
            % time we only count the coupling from the cost)
for i = 1:Nx*tFIR+Nu*(tFIR-1)
    indices{i} = [];
    for j = 1:Nx*tFIR+Nu*(tFIR-1)
        if C(i,j) ~= 0
            indices{i} = [indices{i} j];
        end
        indices{i} = unique(indices{i});
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
    
Psi_prev = ones(Nx*tFIR + Nu*(tFIR-1),Nx); % Just so the while doesn't break

rho = 300;

count = 0; conv = [1];
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

    mu = 50; 

        while norm(conv1)>10^(-3) || norm(conv2)>10^(-4)


        %% Phi & X

        sys = 1;
            for i = r{sys}

                n = max(length(s_r{sys}(find(r{sys}==i),:)));
                m = max(length(indices{i}));

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
                    b{k} = Z_admm{k}-Y{i}{k};
                    Mjb_sum = Mjb_sum + M{j}'*b{k};
                end


                if i <= Nx*tFIR && i >= Nx*2
                    K_proxi = K(2*i-1:2*i,indices{i});
                    model.Q = sparse((C_proxi*M2)'*(C_proxi*M2)+rho/2*(M1'*M1)+mu/2*(Mj_sum));
                    model.A = sparse(K_proxi*M2); 
                    model.obj = -rho*a*M1 -mu*Mjb_sum';
                    model.rhs = [up up]; 
                    model.sense = '<';
                    model.lb = -100*ones(m+n-1,1);
                    params.outputflag = 0;

%                     if t > 1
%                     tic;
%                     end
                    result = gurobi(model, params);
                    W = result.x(:); 
                    if t>1
    %                 [~] = toc;
    %                 time(index) = time(index) + toc;
                      times(index) = times(index) + result.runtime;
                    end

                else
                    W = pinv(2*((C_proxi*M2)'*C_proxi*M2)+rho*(M1'*M1)+mu*Mj_sum)*(rho*M1'*a'+mu*Mjb_sum);
                end

                Phi_loc{i} = (M1*W)';

                X_i = M2*W;
                X{i} = zeros(Nx*tFIR+Nu*(tFIR-1),1); X{i}(indices{i}) = X_i;

            end

        for sys = 2:Nx

            for i = r{sys}

                n = max(length(s_r{sys}(find(r{sys}==i),:)));
                m = max(length(indices{i}));

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
                    b{k} = Z_admm{k}-Y{i}{k};
                    Mjb_sum = Mjb_sum + M{j}'*b{k};
                end


                if i <= Nx*tFIR && i >= Nx*2
                    K_proxi = K(2*i-1:2*i,indices{i});
                    model.Q = sparse((C_proxi*M2)'*(C_proxi*M2)+rho/2*(M1'*M1)+mu/2*(Mj_sum));
                    model.A = sparse(K_proxi*M2); 
                    model.obj = -rho*a*M1 -mu*Mjb_sum';
                    model.rhs = [up up]; 
                    model.sense = '<';
                    model.lb = -100*ones(m+n-1,1);
                    params.outputflag = 0;

                    result = gurobi(model, params);
                    W = result.x(:); 

                else
                    W = pinv(2*((C_proxi*M2)'*C_proxi*M2)+rho*(M1'*M1)+mu*Mj_sum)*(rho*M1'*a'+mu*Mjb_sum);
                end

                Phi_loc{i} = (M1*W)';

                X_i = M2*W;
                X{i} = zeros(Nx*tFIR+Nu*(tFIR-1),1); X{i}(indices{i}) = X_i;

            end
        end


        %% Z

        Z_old = Z_admm;

        sys = 1; % Separate since we want to measure time
        if t > 1
        tic;
        end
        for j = r{sys}
            Z_admm{j} = 0;
            for i = indices{j}
                Z_admm{j} = Z_admm{j} + (X{i}(j)+Y{i}{j})/length(indices{j});
            end
        end
        if t > 1
        [~] = toc;
        times(index) = times(index) + toc;
        end

        for sys = 2:Nx
            for j = r{sys}
                Z_admm{j} = 0;
                for i = indices{j}
                    Z_admm{j} = Z_admm{j} + (X{i}(j)+Y{i}{j})/length(indices{j});
                end
            end
        end

        %% Lagrange multiplier

        sys = 1; % Separate since we want to measure time
        if t > 1
        tic;
        end
        for i = r{sys}
            for j = indices{i}
                Y{i}{j} = Y{i}{j} + X{i}(j) - Z_admm{j};
            end
        end
        if t > 1
        [~] = toc;
        times(index) = times(index) + toc;
        end

        for sys = 2:Nx
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
        if counter > 500
            break
            disp('Consensus ADMM did not converge')
        end

        end
    
        if t > 1
            iters(index) = iters(index) + counter;
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
    
    % Solve for each column
    i = 1; % Separate the first row since we want to measure time
    IZA_ZB_loc = IZA_ZB(:,s_c{i}); row_all_zeros = find(all(IZA_ZB_loc == 0,2)); keep_indices = setdiff(linspace(1,Nx*tFIR,Nx*tFIR),row_all_zeros);
    IZA_ZB_loc = IZA_ZB(keep_indices,s_c{i}); E1_loc = E1(keep_indices,c{i});
    if t > 1
        tic
    end
    AUX_matrix = IZA_ZB_loc'*pinv(IZA_ZB_loc*IZA_ZB_loc');
    Psi_loc{i} = (Phi_loc_col{i}+Lambda_loc_col{i})+AUX_matrix*(E1_loc-IZA_ZB_loc*(Phi_loc_col{i}+Lambda_loc_col{i}));
    if t > 1
        [~] = toc;
        times(index) = times(index) + toc;
    end
    
    for i = 2:Nx
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
    count = count + 1
    if count > 5000
        disp ('ADMM did not converge')
        break
    end

end

    if t > 1
        iters(index) = iters(index) + count;
    end

%% Dynamics
    
    % Compute the control action (in a localized way)
        u(:,t) = Phi(1+Nx*tFIR:Nx*tFIR+Nu,:)*xi;
       
    % Simulate what the dynamics are given that action
    x(:,t+1) = Phi(1+Nx:2*Nx,:)*xi; % Since there is no noise x_ref = x
    
    % Update the initial condition
    xi = x(:,t+1);
    
end
times(index) = times(index)/(tSim-1);

%% Validation (include in loop)
% Centralized MPC (for validation + comparison)
coupling = true;
[xVal, uVal, ~] = mpc_centralized(Nx, Nu, A, B, d, Q, S, tFIR, tSim, x0, up, KSmall, coupling);
   
end

%% Plot
figure(1)
subplot(1,2,1)
plot(localities, times,'m-s','LineWidth',2)
hold on
plot(localities, timeCents,'b-s','LineWidth',2)
xlabel('$$\#\ pendulums\ in\ the\ network$$','Interpreter','latex','Fontsize', 10)
ylabel('$$Avg\ runtime\ per\ MPC\ iteration\ for\ each\ state\ (s)$$','Interpreter','latex','Fontsize', 10)
leg4 = legend('$$Localized\ ADMM\ Solution$$', '$$Centralized\ Solution$$');
set(leg4,'Interpreter','latex','Fontsize', 8);

subplot(1,2,2)
plot(localities, iters,'m-s','LineWidth',2)
xlabel('$$\#\ pendulums\ in\ the\ network$$','Interpreter','latex','Fontsize', 10)
ylabel('$$Avg\ \#\ ADMM\ iters\ per\ MPC\ iteration\ for\ each\ state\ (s)$$','Interpreter','latex','Fontsize', 10)
