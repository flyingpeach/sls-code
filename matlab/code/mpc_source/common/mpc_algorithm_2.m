function [x, u, avgTime, avgIter] = mpc_algorithm_2(sys, x0, params, ...
                                    indices, C)
% Inputs
%   sys     : LTISystem containing system matrices (A, B2) and Nx, Nu 
%   x0      : Initial system state
%   params  : MPCParams containing parameters for mpc
% Outputs
%   x       : State
%   u       : Input
%   avgTime : Average runtime (Steps 4+6+8+11) per MPC iteration for each state
%   avgIter : Average ADMM iters per MPC iteration for each state

% Note that we include the first state as a representative per-state 
% measurement for runtime calculations. Both time and iteration
% calculations exclude the first mpc iteration to omit warmup effects

%% Setup
fprintf('MPC Algorithm 2\n')
params.sanity_check_alg_2();

% For ease of notation
Nx = sys.Nx; Nu = sys.Nu;
locality = params.locality_;
tFIR     = params.tFIR_;
tHorizon = params.tHorizon_;

maxIters = params.maxIters_;
rho      = params.rho_;
eps_d    = params.eps_d_;
eps_p    = params.eps_p_;

maxItersCons = params.maxItersCons_;
mu           = params.mu_;
eps_x        = params.eps_x_;
eps_z        = params.eps_z_;

% ADMM variables
Phi    = zeros(Nx*tFIR + Nu*(tFIR-1),Nx);
Psi    = zeros(Nx*tFIR + Nu*(tFIR-1),Nx);
Lambda = zeros(Nx*tFIR + Nu*(tFIR-1),Nx);
Y_locs = cell(1, Nx*tFIR+Nu*(tFIR-1));
Z_locs = cell(1, Nx*tFIR+Nu*(tFIR-1));

for i = 1:Nx*tFIR+Nu*(tFIR-1) % Initialize Y, Z
    if ~isempty(indices{i})
        Z_locs{i} = 0;
        for j = indices{i}
            Y_locs{i}{j} = 0;
        end
    else
    end
end

% State + control
x       = zeros(Nx, tHorizon);
u       = zeros(Nu, tHorizon);
x(:, 1) = x0;

% Track time & iterations
totalTime = 0;
totalIter = 0;

% SLS constraints
Eye = [eye(Nx); zeros(Nx*(tFIR-1),Nx)];
ZAB = get_sls_constraint(sys, tFIR);

% Locality setup
[r_loc, m_loc] = get_r_m_locality(sys, locality, tFIR);
[c, s_c]       = get_column_locality(sys, tFIR, r_loc, m_loc);
[r, s_r]       = get_row_locality(sys, tFIR, r_loc, m_loc);

%% MPC
for t = 1:tHorizon
    fprintf('Calculating time %d of %d\n', t, tHorizon); % display progress
    x_t = x(:,t); % update initial condition
    
    for iter=1:maxIters % ADMM outer loop
        Psi_prev = Psi;
        
        % Separate Psi, Lambda into rows
        if params.solnMode_ == MPCSolMode.ClosedForm
            [Psi_rows, Lambda_rows] = separate_rows_2(sys, tFIR, r, s_r, r_loc, m_loc, Psi, Lambda);
        elseif params.solnMode_ == MPCSolMode.Explicit
            % TODO
        elseif params.solnMode_ == MPCSolMode.UseSolver
            % TODO
        end
        
        for consIter=1:maxItersCons
            Z_prev_locs = Z_locs;

            % Step 4: Solve (20a) to get local Phi, X            
            Phi_locs = cell(1, Nx*tFIR + Nu*(tFIR-1));
            X_locs   = cell(1, Nx*tFIR + Nu*(tFIR-1));

            if params.solnMode_ == MPCSolMode.ClosedForm
                for i_ = 1:Nx
                    if t > 1 && i_ == 1; tic; end

                    for i = r{i_}
                        n     = max(length((s_r{i_}(find(r{i_}==i),:))));
                        x_ri  = x_t(s_r{i_}(find(r{i_}==i),:));
                        i_new = find(indices{i} == i);
                        ci    = C(i, indices{i})
                        
                        [Phi_locs{i}, X_locs{i}] = eqn_20a_closed(x_ri, Psi_rows{i}, Lambda_rows{i}, ...
                                                                  Z_locs, y_rowi, indices{i}, i_new, ci, n, rho);                        
                    end
                    if t > 1 && i_ == 1; totalTime = totalTime + toc; end
                end
            elseif params.solnMode_ == MPCSolMode.Explicit
                % TODO
            elseif params.solnMode_ == MPCSolMode.UseSolver
                % TODO
            end
            
            % Step 6: Update Z (Step 5 implicitly done in this step)
            for i_ = 1:Nx
                if t > 1 && i_ == 1; tic; end
                for i = r{i_}
                    Z_locs{i} = 0;
                    for j = indices{i}
                        Z_locs{i} = Z_locs{i} + (X_locs{j}(i)+Y_locs{j}{i})/length(indices{i});
                    end
                end 
                if t > 1 && i_ == 1; totalTime = totalTime + toc; end
            end
       
            % Step 8: Update Y (Step 7 implicitly done in this step)            
            for i_ = 1:Nx
                if t > 1 && i_ == 1; tic; end
                for i = r{i_}
                    for j = indices{i}
                        Y_locs{i}{j} = Y_locs{i}{j} + X_locs{i}(j) - Z_locs{j};
                    end
                end
                if t > 1 && i_ == 1; totalTime = totalTime + toc; end
            end
            
            % Step 9: Check convergence of inner loop
            converged = true;
            
            
            % uncleaned code below this line --------------------
            
            
            % Consensus convergence criterium
            primal = zeros(1, Nx);
            for i_ = 1:Nx
                for i = r{i_}
                    average = zeros(Nx*tFIR+Nu*(tFIR-1),1);
                    for j = indices{i}
                        average(j) = Z_locs{j};
                    end
                    primal(i) = norm(X_locs{i}-average, 'fro');
                end
            end
            primalRes = max(primal); % primal residue

            dual = zeros(1, Nx);
            for i_ = 1:Nx
                for j = r{i_}
                    dual(j) = norm(Z_locs{j}-Z_prev_locs{j}, 'fro');
                end
            end
            dualRes = max(dual); % dual residue

            consCriterion_failed = false;
            if primalRes > eps_x || dualRes > eps_z
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
        for i_ = 1:Nx
            for i = r{i_}
                Phi(i,s_r{i_}(find(r{i_}==i),:)) = Phi_locs{i};
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
        for i_ = 1:Nx
              local_phi = Phi(r{i_},s_r{i_}(tFIR,:));
              local_psi = Psi(r{i_},s_r{i_}(tFIR,:));
              local_psi_prev = Psi_prev(r{i_},s_r{i_}(tFIR,:));

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

avgTime = totalTime / (tHorizon - 1);
avgIter = totalIter / (tHorizon - 1);

end