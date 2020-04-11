function [x, u, avgTime, avgIter] = mpc_algorithm1_gurobi(Nx, Nu, A, B, d, tFIR, tSim, x0, ... % system
                              eps_d, eps_p, rho, maxIters, ... % admm
                              up, low) % gurobi

% ADMM variables
Phi    = zeros(Nx*tFIR + Nu*(tFIR-1),Nx);
Psi    = zeros(Nx*tFIR + Nu*(tFIR-1),Nx);
Lambda = zeros(Nx*tFIR + Nu*(tFIR-1),Nx);

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
        
        Phi_loc = cell(1, Nx);
        sys = 1;
        for i = r{sys}
            n = max(length((s_r{sys}(find(r{sys}==i),:))));
            xi_i = x_t(s_r{sys}(find(r{sys}==i),:));
            if i <= Nx*tFIR
                model.Q = sparse(xi_i*xi_i'+rho/2*eye(n));
                model.obj = rho*(-Psi_loc_row{i}+Lambda_loc_row{i});
                model.A = sparse([xi_i'; xi_i']);
                model.rhs = [up low];
                model.sense = '<>';
                model.lb = -100*ones(n,1);
                params.outputflag = 0;
                
                result = gurobi(model, params);
                if t>1
                  totalTime = totalTime + result.runtime;
                end
                Phi_loc{i} = result.x(:);
                
            else
                ADMM_matrix = inv(2*xi_i*xi_i'+rho*eye(n));
                Phi_loc{i} = rho*(Psi_loc_row{i}-Lambda_loc_row{i})*ADMM_matrix;
            end
        end
        
        for sys = 2:Nx
            for i = r{sys}
                n = max(length((s_r{sys}(find(r{sys}==i),:))));
                xi_i = x_t(s_r{sys}(find(r{sys}==i),:));
                if i <= Nx*tFIR
                    model.Q = sparse(xi_i*xi_i'+rho/2*eye(n));
                    model.obj = rho*(-Psi_loc_row{i}+Lambda_loc_row{i});
                    model.A = sparse([xi_i'; xi_i']);
                    model.rhs = [up low];
                    model.sense = '<>';
                    model.lb = -100*ones(n,1);
                    params.outputflag = 0;
                    
                    result = gurobi(model, params);
                    Phi_loc{i} = result.x(:);
                    
                else
                    ADMM_matrix = inv(2*xi_i*xi_i'+rho*eye(n));
                    Phi_loc{i} = rho*(Psi_loc_row{i}-Lambda_loc_row{i})*ADMM_matrix;
                end
            end
        end
        
        % Build Phi matrix
        for sys = 1:Nx
            for i = r{sys}
                Phi(i,s_r{sys}(find(r{sys}==i),:)) = Phi_loc{i};
            end
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
        % for timing, separate first column
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
         
        % Build Phi matrix
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
            break; % convergence criterion passed, ex_tt admm iterations
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