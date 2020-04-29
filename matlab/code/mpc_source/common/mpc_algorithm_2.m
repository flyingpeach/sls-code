function [x, u, avgTime, avgIter] = mpc_algorithm_2(sys, x0, params)
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

up = params.constrUpperbnd_;

% ADMM variables
Phi    = zeros(Nx*tFIR + Nu*(tFIR-1),Nx);
Psi    = zeros(Nx*tFIR + Nu*(tFIR-1),Nx);
Lambda = zeros(Nx*tFIR + Nu*(tFIR-1),Nx);
Y_locs = cell(1, Nx*tFIR+Nu*(tFIR-1));
Z_locs = cell(1, Nx*tFIR+Nu*(tFIR-1));

Ksmall  = params.constrMtx_;
% Build the big constraint matrix
K = [zeros(size(Ksmall)) zeros(2*Nx,(tFIR-1)*Nx)];
K = [K; zeros(2*Nx,Nx) zeros(size(Ksmall)) zeros(2*Nx,(tFIR-2)*Nx)];
for t = 2:tFIR-1
    K = [K; zeros(2*Nx,t*Nx) Ksmall zeros(2*Nx,(tFIR-t-1)*Nx)];
end

% Build cost matrix (block diagonal)
C = [];
for t = 0:tFIR-1
    C = blkdiag(C, params.Q_);
end
for t = 0:tFIR-2
    C = blkdiag(C, params.R_);
end    

% Note: only coupling from cost matrix is considered
%       coupling from constraints is ignored
indices = cell(1, length(C));
for i = 1:length(C)
    for j = 1:length(C)
        if C(i,j) ~= 0
            indices{i} = [indices{i} j];
        end
    end
end

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
    
    for iter=1:maxIters % ADMM (outer loop)
        Psi_prev = Psi;
        
        % Separate Psi, Lambda into rows
        [Psi_rows, Lambda_rows] = separate_rows_2(sys, tFIR, r, s_r, r_loc, m_loc, Psi, Lambda);
        
        for consIter=1:maxItersCons % ADMM consensus (inner loop)
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
                        ci    = C(i, indices{i});
                        
                        [Phi_locs{i}, X_locs{i}] = eqn_20a_closed(x_ri, Psi_rows{i}, Lambda_rows{i}, ...
                                                                  Z_locs, Y_locs{i}, indices{i}, i_new, ci, n, rho, mu);                        
                    end
                    if t > 1 && i_ == 1; totalTime = totalTime + toc; end
                end
            elseif params.solnMode_ == MPCSolMode.UseSolver
                for i_ = 1:Nx
                    if t > 1 && i_ == 1; tic; end

                    for i = r{i_}
                        n     = max(length((s_r{i_}(find(r{i_}==i),:))));
                        x_ri  = x_t(s_r{i_}(find(r{i_}==i),:));
                        i_new = find(indices{i} == i);
                        ci    = C(i, indices{i});

                        if i <= Nx*tFIR && i >= Nx*2
                            ki = K(2*i-1:2*i, indices{i});
                            [Phi_locs{i}, X_locs{i}, time] = eqn_20a_solver(x_ri, Psi_rows{i}, Lambda_rows{i}, ...
                                                                            Z_locs, Y_locs{i}, ki, indices{i}, i_new, ci, n, rho, mu, up);
                        else
                            [Phi_locs{i}, X_locs{i}] = eqn_20a_closed(x_ri, Psi_rows{i}, Lambda_rows{i}, ...
                                                                      Z_locs, Y_locs{i}, indices{i}, i_new, ci, n, rho, mu);
                        end
                    end
                    if t > 1 && i_ == 1; totalTime = totalTime + time; end
                end                
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
            
            % Step 9: Check convergence of ADMM consensus
            converged = true;           
            for i_ = 1:Nx
                for i = r{i_}
                    z_av = zeros(Nx*tFIR+Nu*(tFIR-1),1);
                    for j = indices{i}
                        z_av(j) = Z_locs{j};
                    end
                    
                    if ~check_convergence_cons(z_av, X_locs{i}, Z_locs{i}, Z_prev_locs{i}, eps_x, eps_z)
                        converged = false;
                        break; % if one fails, can stop checking the rest
                    end
                end
            end
            
            if converged
                break; % exit ADMM consensus iterations
            end
        end
        
        if t > 1; totalIter = totalIter + consIter; end
                
        if ~converged    
            fprintf('ADMM consensus reached %d iters and did not converge\n', maxItersCons);
        end

        % Step 10: Build entire Phi matrix        
        for i_ = 1:Nx
            for i = r{i_}
                Phi(i,s_r{i_}(find(r{i_}==i),:)) = Phi_locs{i};
            end
        end

        % Separate Phi, Lambda into columns
        [Phi_cols, Lambda_cols] = separate_cols(sys, c, s_c, Phi, Lambda);
        
        % Step 11: Solve (16b) to get local Psi
        Psi_locs = cell(1, Nx);
        for i = 1:Nx
            if t > 1 && i == 1; tic; end

            ZABi     = ZAB(:, s_c{i});
            zeroRow  = find(all(ZABi == 0, 2));
            keepIdxs = setdiff(linspace(1,Nx*tFIR,Nx*tFIR), zeroRow);
            ZABi     = ZAB(keepIdxs, s_c{i}); 
            Eyei     = Eye(keepIdxs, c{i});

            Psi_locs{i} = eqn_16b(Phi_cols{i}, Lambda_cols{i}, ZABi, Eyei);

            if t > 1 && i == 1; totalTime = totalTime + toc; end            
        end
        
        % Step 12: Build entire Psi matrix
        for i = 1:Nx
            Psi(s_c{i},c{i}) = Psi_locs{i};
        end

        % Step 13: Update Lambda
        Lambda = Lambda + Phi - Psi; 
            
        % Step 14: Check convergence of ADMM (outer loop)
        converged = true;
        for i = 1:Nx
              phi_loc      = Phi(r{i},s_r{i}(tFIR,:));
              psi_loc      = Psi(r{i},s_r{i}(tFIR,:));
              psi_prev_loc = Psi_prev(r{i},s_r{i}(tFIR,:));

              if ~check_convergence(phi_loc, psi_loc, psi_prev_loc, eps_p, eps_d)
                  converged = false;
                  break; % if one fails, can stop checking the rest
              end
        end
        
        if converged
            break; % exit ADMM iterations
        end
    end

    if t > 1; totalIter = totalIter + iter; end
    
    if ~converged
        fprintf('ADMM reached %d iters and did not converge\n', maxIters);
    end
    
    % Compute control + state
    u(:,t)   = Phi(1+Nx*tFIR:Nx*tFIR+Nu,:)*x_t;
    x(:,t+1) = Phi(1+Nx:2*Nx,:)*x_t; % since no noise, x_ref = x
end

avgTime = totalTime / (tHorizon - 1);
avgIter = totalIter / (tHorizon - 1);
end