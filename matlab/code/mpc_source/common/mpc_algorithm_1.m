function [x, u, avgTime, avgIter] = mpc_algorithm_1(sys, x0, params)
% Inputs
%   sys     : LTISystem containing system matrices (A, B2) and Nx, Nu 
%   x0      : Initial system state
%   params  : MPCParams containing parameters for mpc
% Outputs
%   x       : State
%   u       : Input
%   avgTime : Average runtime (Steps 4+6) per MPC iteration for each state
%   avgIter : Average ADMM iters per MPC iteration for each state

% Note that we include the first state as a representative per-state 
% measurement for runtime calculations. Both time and iteration
% calculations exclude the first mpc iteration to omit warmup effects

%% Setup
fprintf('MPC Algorithm 1\n')
params.sanity_check_alg_1();

% For ease of notation
Nx = sys.Nx; Nu = sys.Nu;
tFIR     = params.tFIR_;
tHorizon = params.tHorizon_;
maxIters = params.maxIters_;
rho      = params.rho_;
eps_d    = params.eps_d_;
eps_p    = params.eps_p_;

% ADMM variables
Phi    = zeros(Nx*tFIR + Nu*(tFIR-1),Nx);
Psi    = zeros(Nx*tFIR + Nu*(tFIR-1),Nx);
Lambda = zeros(Nx*tFIR + Nu*(tFIR-1),Nx);

% State + control
x       = zeros(Nx, tHorizon);
u       = zeros(Nu, tHorizon);
x(:, 1) = x0;

% Track time & iterations
totalTime = 0;
totalIter = 0;

% Other setup
Eye = [eye(Nx); zeros(Nx*(tFIR-1),Nx)];
ZAB = get_sls_constraint(sys, tFIR);

[r_loc, m_loc] = get_r_m_locality(sys, tFIR);
[c, s_c]       = get_column_locality(r_loc, m_loc, tFIR);
[r, s_r]       = get_row_locality(r_loc, m_loc, tFIR);

%% MPC
for t = 1:tHorizon
    fprintf('Calculating time %d of %d\n', t, tHorizon); % display progress
    x_t = x(:,t); % update initial condition

    for iter=1:maxIters % ADMM loop
        Psi_prev = Psi;
        
        [Psi_rows, Lambda_rows] = separate_rows(r, s_r, r_loc, m_loc, tFIR, Nx);
       
        % Step 4: Solve (16a) to get local Phi
        Phi_locs = cell(1, Nx);
        for i = 1:Nx
            if t > 1 && i == 1; tic; end
            
            x_ri        = x_t(s_r{i}(tFIR, :)); % local state
            Phi_locs{i} = eqn_16a_closed(s_r{i}, x_ri, Psi_rows{i}, Lambda_rows{i}, rho);
            
            if t > 1 && i == 1; totalTime = totalTime + toc; end
        end
        
        % Step 5: Build entire Phi matrix
        for i = 1:Nx
            Phi(r{i},s_r{i}(tFIR,:)) = Phi_locs{i};
        end

        [Phi_cols, Lambda_cols] = separate_cols(c, s_c, Nx);
        
        % Step 6: Solve (16b) to get local Psi
        Psi_locs = cell(1, Nx);
        for i = 1:Nx
            if t > 1 && i == 1; tic; end

            Psi_locs{i} = eqn_16b(c{i}, s_c{i}, Phi_cols{i}, Lambda_cols{i}, ZAB, Eye);

            if t > 1 && i == 1; totalTime = totalTime + toc; end            
        end
        
        % Step 7: Build entire Psi matrix
        for i = 1:Nx
            Psi(s_c{i},c{i}) = Psi_locs{i};
        end
                     
        % Step 8: Update Lambda
        Lambda = Lambda + Phi - Psi;
        
        % Step 9: Check convergence
        converged = true;
        for i = 1:Nx
              phi_loc      = Phi(r{i},s_r{i}(tFIR,:));
              psi_loc      = Psi(r{i},s_r{i}(tFIR,:));
              psi_prev_loc = Psi_prev(r{i},s_r{i}(tFIR,:));

              if ~check_convergence(phi_loc, psi_loc, psi_prev_loc, eps_p, eps_d);
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
    u(:,t) = Phi(1+Nx*tFIR:Nx*tFIR+Nu,:)*x_t;
    x(:,t+1) = Phi(1+Nx:2*Nx,:)*x_t; % since no noise, x_ref = x
end

avgTime = totalTime / (tHorizon - 1);
avgIter = totalIter / (tHorizon - 1);
end