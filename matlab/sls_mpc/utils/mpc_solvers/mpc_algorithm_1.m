function [x, u, time, iters] = mpc_algorithm_1(sys, x0, params)
% Inputs
%   sys     : LTISystem containing system matrices (A, B2) and Nx, Nu 
%   x0      : Starting system state
%   params  : MPCParams containing parameters for mpc
% Outputs
%   x       : Next state if MPC input is used
%   u       : Current input as calculated by MPC
%   time    : Total runtime (Steps 4+6) per state
%   iters   : Total ADMM iters per state
%
% Note that we include the first state as a representative per-state 
% measurement for runtime calculations

%% Setup
params.sanity_check_alg_1();

% For ease of notation
Nx = sys.Nx; Nu = sys.Nu;
locality = params.locality_;
tFIR     = params.tFIR_;
maxIters = params.maxIters_;
rho      = params.rho_;

nVals  = Nx*tFIR + Nu*(tFIR-1);
% ADMM variables
Phi    = zeros(nVals, Nx);
Psi    = zeros(nVals, Nx);
Lambda = zeros(nVals, Nx);

% Cost matrices (both are diagonal)
C = build_cost_mtx(params);
K = build_constr_mtx(sys, params);

% Track runtime
time = 0;

% SLS constraints
Eye = [eye(Nx); zeros(Nx*(tFIR-1),Nx)];
ZAB = get_sls_constraint(sys, tFIR);

% Get indices corresponding to rows / columns / localities
[r_loc, m_loc] = get_r_m_locality(sys, locality);
[c, s_c]       = get_col_locality(sys, tFIR, r_loc, m_loc);
[r, s_r]       = get_row_locality(sys, tFIR, r_loc, m_loc);

% Column-wise partition SLS constraints and pre-calculate inverse
zabs  = cell(Nx, 1);
eyes  = cell(Nx, 1);
zabis = cell(Nx, 1); % inverse

for i = 1:Nx
    zab_     = ZAB(:, s_c{i});
    zeroRows = find(all(zab_ == 0, 2));
    keepRows = setdiff(1:tFIR*Nx, zeroRows);           
    zabs{i}  = ZAB(keepRows, s_c{i}); 
    eyes{i}  = Eye(keepRows, c{i});
    
    zabis{i} = zabs{i}'*pinv(zabs{i}*zabs{i}');
end

%% MPC
for iters=1:maxIters % ADMM loop
    Psi_prev = Psi;
        
    % Separate Psi, Lambda into rows (with sparsity)
    Psi_rows    = separate_rows(sys, tFIR, r, s_r, Psi);
    Lambda_rows = separate_rows(sys, tFIR, r, s_r, Lambda);
      
    % Step 4: Solve (16a) to get local rows of Phi
    Phi_rows = cell(nVals, 1);
        
    for i = 1:Nx
        if i == 1; tic; end
         
        for j = 1:length(r{i})
            % Note: Alg1 has no coupling, so C, K are diagonal
            row   = r{i}(j);
            x_loc = x0(s_r{i}{j}); % observe local state
            cost_ = C(row, row);
                
            if K(row, row) % has constraint
                if row <= tFIR*Nx % state constraint
                    b1_ = params.stateUB_(i) / K(row, row);
                    b2_ = params.stateLB_(i) / K(row, row);
                else % input constraint
                    inputIdx = find(sys.B2(i,:));
                    b1_ = params.inputUB_(inputIdx) / K(row, row);
                    b2_ = params.inputLB_(inputIdx) / K(row, row);
                end
                b1  = max(b1_,b2_); b2 = min(b1_,b2_); % in case of negative signs
                Phi_rows{row} = eqn_16a_explicit(x_loc, Psi_rows{row}, Lambda_rows{row}, b1, b2, cost_, rho);
            else % no constraint, use closed form
                Phi_rows{row} = eqn_16a_closed(x_loc, Psi_rows{row}, Lambda_rows{row}, cost_, rho); 
            end
        end
            
        if i == 1; time = time + toc; end 
    end
        
    % Step 5: Build entire Phi matrix
    Phi = build_from_rows(sys, r, s_r, Phi_rows, size(Phi));

    % Separate Phi, Lambda into columns
    Phi_cols    = separate_cols(sys, c, s_c, Phi);
    Lambda_cols = separate_cols(sys, c, s_c, Lambda);
        
    % Step 6: Solve (16b) to get local columns of Psi
    Psi_cols = cell(Nx, 1);
    for i = 1:Nx
        if i == 1; tic; end
        Psi_cols{i} = eqn_16b(Phi_cols{i}, Lambda_cols{i}, zabs{i}, eyes{i}, zabis{i});
        if i == 1; time = time + toc; end            
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
    
if ~converged
    fprintf('ADMM reached %d iters and did not converge\n', maxIters);
end
    
% Compute control + state
u = Phi(1+Nx*tFIR:Nx*tFIR+Nu,:)*x0;
x = Phi(1+Nx:2*Nx,:)*x0; % since no noise, x_ref = x

end