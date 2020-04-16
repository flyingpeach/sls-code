function [x, u, avgTime, avgIter] = mpc_algorithm2(Nx, Nu, A, B, d, tFIR, tSim, x0, ... % system
                              eps_d, eps_p, rho, maxIters, ...
                              indices, eps_pres, eps_dres, mu, ...
                              maxConsensusIters, C)
% ADMM variables
Phi = zeros(Nx*tFIR + Nu*(tFIR-1),Nx);
Psi = zeros(Nx*tFIR + Nu*(tFIR-1),Nx);
Lambda = zeros(Nx*tFIR + Nu*(tFIR-1),Nx);

Z_admm = cell(1, Nx*tFIR+Nu*(tFIR-1));
Y      = cell(1, Nx*tFIR+Nu*(tFIR-1));

for i = 1:Nx*tFIR+Nu*(tFIR-1)
    if ~isempty(indices{i})
        Z_admm{i} = 0;
        for j = indices{i}
            Y{i}{j} = 0;
        end
    else
    end
end

% Keep track of state + control
x       = zeros(Nx, tSim);
u       = zeros(Nu, tSim);
x(:, 1) = x0;

% Set up constraints
[E1, IZA_ZB] = get_sls_constraints(tFIR, Nx, Nu, A, B);
[r, c, s_r, s_c, LocalityR, LocalityM] = get_locality_constraints(tFIR, Nx, Nu, A, B, d);

% Track time / iterations
totalTime = 0;
totalIter = 0;

for t = 1:tSim
    fprintf('Calculating time %d of %d\n', t, tSim); % display progress
    x_t = x(:,t); % update initial condition
    
    for iter=1:maxIters
        Psi_prev = Psi;
        
        Psi_loc_row    = cell(1, Nx);
        Lambda_loc_row = cell(1, Nx);
        
        % Separate into rows
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
        
        % Step 4 of Algorithm 2
        for consIter=1:maxConsensusIters
            
            Phi_loc = cell(1, Nx*tFIR + Nu*(tFIR-1));
            X       = cell(Nx*tFIR + Nu*(tFIR-1));
            % Update Phi and X    
            sys = 1;
            tic;
            for i = r{sys}
                n = max(length(s_r{sys}(find(r{sys}==i),:)));

                % Define parameters
                m_j = max(length(indices{i}));
                C_proxi = C(i,indices{i});
                M1 = [eye(n) zeros(n,m_j-1)];
                Id = eye(m_j-1);
                xi_i = x_t(s_r{sys}(find(r{sys}==i),:));

                i_new = find(indices{i} == i);

                M2 = [zeros(i_new-1,n) Id(1:i_new-1,:); xi_i' zeros(1,m_j-1); zeros(m_j-i_new,n) Id(i_new:end,:)];
                a = Psi_loc_row{i}-Lambda_loc_row{i};

                Mj_sum = 0; Mjb_sum = 0;
                   
                M = cell(1, m_j);
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
            totalTime = totalTime + toc;
            
            for sys = 2:Nx
                for i = r{sys}

                    n = max(length(s_r{sys}(find(r{sys}==i),:)));

                    % Define parameters
                    m_j = max(length(indices{i}));
                    C_proxi = C(i,indices{i});
                    M1 = [eye(n) zeros(n,m_j-1)];
                    Id = eye(m_j-1);
                    xi_i = x_t(s_r{sys}(find(r{sys}==i),:));

                    i_new = find(indices{i} == i);

                    M2 = [zeros(i_new-1,n) Id(1:i_new-1,:); xi_i' zeros(1,m_j-1); zeros(m_j-i_new,n) Id(i_new:end,:)];
                    a = Psi_loc_row{i}-Lambda_loc_row{i};

                    Mj_sum = 0; Mjb_sum = 0;
                    
                    M = cell(1, m_j);
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
            
            % Update Z
            Z_old = Z_admm;
            
            sys = 1;
            tic;
            for j = r{sys}
                Z_admm{j} = 0;
                for i = indices{j}
                    Z_admm{j} = Z_admm{j} + (X{i}(j)+Y{i}{j})/length(indices{j});
                end
            end 
            totalTime = totalTime + toc;            
            
            for sys = 2:Nx
                for j = r{sys}
                    Z_admm{j} = 0;
                    for i = indices{j}
                        Z_admm{j} = Z_admm{j} + (X{i}(j)+Y{i}{j})/length(indices{j});
                    end
                end 
            end
        
            % Lagrange multiplier
            sys = 1;
            tic;
            for i = r{sys}
                for j = indices{i}
                    Y{i}{j} = Y{i}{j} + X{i}(j) - Z_admm{j};
                end
            end
            totalTime = totalTime + toc;
            
            for sys = 2:Nx
                for i = r{sys}
                    for j = indices{i}
                        Y{i}{j} = Y{i}{j} + X{i}(j) - Z_admm{j};
                    end
                end
            end
            
            % Consensus convergence criterium
            primal = zeros(1, Nx);
            for sys = 1:Nx
                for i = r{sys}
                    average = zeros(Nx*tFIR+Nu*(tFIR-1),1);
                    for j = indices{i}
                        average(j) = Z_admm{j};
                    end
                    primal(i) = norm(X{i}-average, 'fro');
                end
            end
            primalRes = max(primal); % primal residue

            dual = zeros(1, Nx);
            for sys = 1:Nx
                for j = r{sys}
                    dual(j) = norm(Z_admm{j}-Z_old{j}, 'fro');
                end
            end
            dualRes = max(dual); % dual residue

            consCriterion_failed = false;
            if primalRes > eps_pres || dualRes > eps_dres
                consCriterion_failed = true;
            end
            
            if ~consCriterion_failed
                break;
            end
        end
                
        if consCriterion_failed        
            fprintf('ADMM consensus reached %d iters and did not converge\n', maxConsensusIters);
        end

        if t > 1
            totalIter = totalIter + consIter;
        end
        
        % Build the big matrix
        for sys = 1:Nx
            for i = r{sys}
                Phi(i,s_r{sys}(find(r{sys}==i),:)) = Phi_loc{i};
            end
        end

        % Separate into columns
        Phi_loc_col = cell(1, Nx);
        Lambda_loc_col = cell(1, Nx);
        for i = 1:Nx
            Phi_loc_col{i} = Phi(s_c{i},c{i});
            Lambda_loc_col{i} = Lambda(s_c{i},c{i});
        end

        Psi_loc = cell(1, Nx);
        
        i = 1;
        IZA_ZB_loc = IZA_ZB(:,s_c{i}); row_all_zeros = find(all(IZA_ZB_loc == 0,2)); keep_indices = setdiff(linspace(1,Nx*tFIR,Nx*tFIR),row_all_zeros);
        IZA_ZB_loc = IZA_ZB(keep_indices,s_c{i}); E1_loc = E1(keep_indices,c{i}); 
        tic;
        AUX_matrix = IZA_ZB_loc'*pinv(IZA_ZB_loc*IZA_ZB_loc');
        Psi_loc{i} = (Phi_loc_col{i}+Lambda_loc_col{i})+AUX_matrix*(E1_loc-IZA_ZB_loc*(Phi_loc_col{i}+Lambda_loc_col{i}));
        totalTime = totalTime + toc;
        
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

        % Lagrange multiplier
        Lambda = Lambda + Phi - Psi; 
            
        % Check convergence locally 
        criterion_failed = false;
        for sys = 1:Nx
              local_phi = Phi(r{sys},s_r{sys}(tFIR,:));
              local_psi = Psi(r{sys},s_r{sys}(tFIR,:));
              local_psi_prev = Psi_prev(r{sys},s_r{sys}(tFIR,:));

              local_diff_d = norm(local_psi-local_psi_prev,'fro');
              local_diff_p = norm(local_phi-local_psi,'fro');

              if local_diff_p > eps_p || local_diff_d > eps_d
                  criterion_failed = true;
                  break; % if one fails, can stop checking the rest
              end
        end

        if ~criterion_failed
            break; % convergence criterion passed, exit admm iterations
        end
    end

    if t > 1
        totalIter = totalIter + iter;
    end
    
    if criterion_failed
        fprintf('ADMM reached %d iters and did not converge\n', maxIters);
    end
    
    % Compute control + state
    u(:,t) = Phi(1+Nx*tFIR:Nx*tFIR+Nu,:)*x_t;
    x(:,t+1) = Phi(1+Nx:2*Nx,:)*x_t; % since no noise, x_ref = x
end

avgTime = totalTime / (tSim - 1);
avgIter = totalIter / (tSim - 1);

end