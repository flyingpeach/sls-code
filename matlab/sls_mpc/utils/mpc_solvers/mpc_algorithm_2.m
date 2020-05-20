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
fprintf('Distributed MPC, Algorithm 2\n')
params.sanity_check_alg_2();

% For ease of notation
Nx = sys.Nx; Nu = sys.Nu;
locality = params.locality_;
tFIR     = params.tFIR_;
tHorizon = params.tHorizon_;

maxIters = params.maxIters_;
rho      = params.rho_;

maxItersCons = params.maxItersCons_;
mu           = params.mu_;

nVals = Nx*tFIR + Nu*(tFIR-1);
% ADMM variables
Phi    = zeros(nVals, Nx);
Psi    = zeros(nVals, Nx);
Lambda = zeros(nVals, Nx);
Y_rows = cell(nVals, 1);
Z_rows = cell(nVals, 1);

% Constraints, costs, and coupling info

% TODO: assumes no input constraints
if params.has_state_cons()
    K = build_constr_mtx(sys, params);
end

C     = build_cost_mtx(params);
cpIdx = get_coupling_indices(C);

for row = 1:nVals % Initialize Y, Z
    if ~isempty(cpIdx{row})
        Z_rows{row} = 0;
        for k = cpIdx{row}
            Y_rows{row}{k} = 0;
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

% Get indices corresponding to rows / columns / localities
[r_loc, m_loc] = get_r_m_locality(sys, locality);
[c, s_c]       = get_col_locality(sys, tFIR, r_loc, m_loc);
[r, s_r]       = get_row_locality(sys, tFIR, r_loc, m_loc);

%% MPC
for t = 1:tHorizon
    fprintf('Calculating time %d of %d\n', t, tHorizon); % display progress
    x_t = x(:,t); % update initial condition
    
    for iter=1:maxIters % ADMM (outer loop)
        Psi_prev = Psi;
        
        % Separate Psi, Lambda into rows (with sparsity)
        Psi_rows    = separate_rows(sys, tFIR, r, s_r, Psi);
        Lambda_rows = separate_rows(sys, tFIR, r, s_r, Lambda);
        
        for consIter=1:maxItersCons % ADMM consensus (inner loop)
            Z_prev_rows = Z_rows;

            % Step 4: Solve (20a) to get local Phi, X            
            Phi_rows = cell(nVals, 1);
            X_rows   = cell(nVals, 1);
            
            for i = 1:Nx
                if t > 1 && i == 1; tic; end
                
                for j = 1:length(r{i})
                    row   = r{i}(j);
                    x_loc = x_t(s_r{i}{j});          % observe local state
                    self_ = find(cpIdx{row} == row); % index of "self-coupling" term
                    cost_ = C(row, cpIdx{row});      % subset of cost matrix
                    
                    isState   = row <= tFIR*Nx; % row represents state
                    % TODO: actually should check if constraint matrix K is
                    % empty
                    stateCons = isState && params.has_state_cons() && ~isempty(params.stateConsMtx_(i,:));

                    % TODO: assumes no input cons
                    if stateCons
                        k_ = K(2*row-1:2*row, cpIdx{row});
                        ub = params.stateUB_;
                        lb = params.stateLB_;
                        
                        [Phi_rows{row}, X_rows{row}] = eqn_20a_solver(x_loc, Psi_rows{row}, Lambda_rows{row}, Y_rows{row}, Z_rows, ...
                                                                      k_, cost_, cpIdx{row}, self_, params, ub, lb);
                    else
                        [Phi_rows{row}, X_rows{row}] = eqn_20a_closed(x_loc, Psi_rows{row}, Lambda_rows{row}, Y_rows{row}, Z_rows, ...
                                                                      cost_, cpIdx{row}, self_, params);
                    end
                end
                
                if t > 1 && i == 1; totalTime = totalTime + toc; end             
            end
               
            % Step 6: Update Z (Step 5 implicitly done in this step)
            for i = 1:Nx
                if t > 1 && i == 1; tic; end
                for row = r{i}
                    Z_rows{row} = 0;
                    for k = cpIdx{row}
                        Z_rows{row} = Z_rows{row} + (X_rows{k}(row)+Y_rows{k}{row})/length(cpIdx{row});
                    end
                end
                if t > 1 && i == 1; totalTime = totalTime + toc; end
            end
       
            % Step 8: Update Y (Step 7 implicitly done in this step)            
            for i = 1:Nx
                if t > 1 && i == 1; tic; end
                for row = r{i}
                    for k = cpIdx{row}
                        Y_rows{row}{k} = Y_rows{row}{k} + X_rows{row}(k) - Z_rows{k};
                    end
                end
                if t > 1 && i == 1; totalTime = totalTime + toc; end
            end
            
            % Step 9: Check convergence of ADMM consensus
            converged = true;           
            for i = 1:Nx
                for row = r{i}
                    z_cp = zeros(nVals, 1);
                    for k = cpIdx{row}
                        z_cp(k) = Z_rows{k};
                    end
                    
                    if ~check_convergence_cons(z_cp, X_rows{row}, Z_rows{row}, Z_prev_rows{row}, params)
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
        Phi = build_from_rows(sys, r, s_r, Phi_rows, size(Phi));

        % Separate Phi, Lambda into columns
        Phi_cols    = separate_cols(sys, c, s_c, Phi);
        Lambda_cols = separate_cols(sys, c, s_c, Lambda);
                
        % Step 11: Solve (16b) to get local Psi
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

        % Step 12: Build entire Psi matrix
        for i = 1:Nx
            Psi(s_c{i}, c{i}) = Psi_cols{i};
        end

        % Step 13: Update Lambda
        Lambda = Lambda + Phi - Psi; 
            
        % Step 14: Check convergence of ADMM (outer loop)
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