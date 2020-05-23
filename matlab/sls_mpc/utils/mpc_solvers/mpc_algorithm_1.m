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
fprintf('Distributed MPC, Algorithm 1\n')
params.sanity_check_alg_1();

% For ease of notation
Nx = sys.Nx; Nu = sys.Nu;
locality = params.locality_;
tFIR     = params.tFIR_;
tHorizon = params.tHorizon_;
maxIters = params.maxIters_;
rho      = params.rho_;

nVals  = Nx*tFIR + Nu*(tFIR-1);
% ADMM variables
Phi    = zeros(nVals, Nx);
Psi    = zeros(nVals, Nx);
Lambda = zeros(nVals, Nx);

% Cost matrix
C = build_cost_mtx(params);

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

% Get indices corresponding to rows / columns / localities
[r_loc, m_loc] = get_r_m_locality(sys, locality);
[c, s_c]       = get_col_locality(sys, tFIR, r_loc, m_loc);
[r, s_r]       = get_row_locality(sys, tFIR, r_loc, m_loc);

%% MPC
for t = 1:tHorizon
    fprintf('Calculating time %d of %d\n', t, tHorizon); % display progress
    x_t = x(:,t); % update initial condition

    for iter=1:maxIters % ADMM loop
        Psi_prev = Psi;
        
        % Separate Psi, Lambda into rows (with sparsity)
        Psi_rows    = separate_rows(sys, tFIR, r, s_r, Psi);
        Lambda_rows = separate_rows(sys, tFIR, r, s_r, Lambda);
      
        % Step 4: Solve (16a) to get local rows of Phi
        Phi_rows = cell(nVals, 1);
        
        for i = 1:Nx
            if t > 1 && i == 1; tic; end
            
            for j = 1:length(r{i})
                row   = r{i}(j);
                x_loc = x_t(s_r{i}{j}); % observe local state
                
                % no coupling in algorithm 1; so C is diagonal
                cost_ = C(row, row);
                
                isState   = row <= tFIR*Nx; % row represents state
                stateCons = isState && params.has_state_cons() && params.stateConsMtx_(i,i);
 
                % TODO: assumes no input cons
                if stateCons
                    m   = params.stateConsMtx_(i,i);
                    b1_ = params.stateUB_ / m;
                    b2_ = params.stateLB_ / m;
                    b1  = max(b1_,b2_); b2 = min(b1_,b2_); % in case of negative signs
                    Phi_rows{row} = eqn_16a_explicit(x_loc, Psi_rows{row}, Lambda_rows{row}, b1, b2, rho);
                else
                    Phi_rows{row} = eqn_16a_closed(x_loc, Psi_rows{row}, Lambda_rows{row}, cost_, rho); 
                end
            end
            
            if t > 1 && i == 1; totalTime = totalTime + toc; end 
        end
        
        % Step 5: Build entire Phi matrix
        Phi = build_from_rows(sys, r, s_r, Phi_rows, size(Phi));

        % Separate Phi, Lambda into columns
        Phi_cols    = separate_cols(sys, c, s_c, Phi);
        Lambda_cols = separate_cols(sys, c, s_c, Lambda);
        
        % Step 6: Solve (16b) to get local columns of Psi
        Psi_cols = cell(Nx, 1);
        for i = 1:Nx
            if t > 1 && i == 1; tic; end

            % Reduce computation by eliminating zero rows
            zab_     = ZAB(:, s_c{i});
            zeroRows = find(all(zab_ == 0, 2));
            keepRows = setdiff(1:tFIR*Nx, zeroRows);           
            zab_     = ZAB(keepRows, s_c{i}); 
            eye_     = Eye(keepRows, c{i});

            Psi_cols{i} = eqn_16b(Phi_cols{i}, Lambda_cols{i}, zab_, eye_);

            if t > 1 && i == 1; totalTime = totalTime + toc; end            
        end
                
        % Step 7: Build entire Psi matrix
        for i = 1:Nx
            Psi(s_c{i}, c{i}) = Psi_cols{i};
        end
                     
        % Step 8: Update Lambda
        Lambda = Lambda + Phi - Psi;
        
        % Step 9: Check convergence
        converged = true;
        for i = 1:Nx
            phi_      = [];
            psi_      = [];
            psi_prev_ = [];
            for j = 1:length(r{i})
                % Due to dimensionality issues, not stacking rows
                % Instead, just make one huge row
                % (since we're checking Frob norm, doesn't matter)
                phi_      = [phi_, Phi(r{i}(j), s_r{i}{j})];
                psi_      = [psi_, Psi(r{i}(j), s_r{i}{j})];
                psi_prev_ = [psi_prev_, Psi_prev(r{i}(j), s_r{i}{j})];
            end
            
            if ~check_convergence(phi_, psi_, psi_prev_, params)
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