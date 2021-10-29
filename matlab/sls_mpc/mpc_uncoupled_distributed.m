function [x, u, stats, warmStartOut] = mpc_uncoupled_distributed(sys, x0, params, warmStartIn)
% Inputs
%   sys     : LTISystem containing system matrices (A, B2) and Nx, Nu
%   x0      : Initial system state
%   params  : MPCParams containing parameters for mpc
%   (optional): warmStartIn: MPCWarmStart containing Psi, Lambda from
%               previous timestep
% Outputs
%   x        : Next state if MPC input is used
%   u        : Current input as calculated by MPC
%   stats    : MPCStats containing runtime/iteration stats
%              time: Runtime per subsystem; time taken for row update, col
%                    update, and consensus updates if applicable
%              iters: ADMM iters per state (outer loop)
%              consIters: Average ADMM consensus iters per state, per iter
%                         (inner loop). Example: if outer loop ran 3 iters
%                         and inner loop ran 2, 5, and 6 iters, consIters
%                         will be 4.3
%   warmStartOut: MPCWarmStart containing Psi, Lambda to be used for
%                 next time iteration, if applicable

%% Setup
sanity_check_actuation(sys);
sanity_check_coupling(sys, params);
params.sanity_check_dist();

% For ease of notation
Nx = sys.Nx; Nu = sys.Nu; T = params.tFIR_;

% Runtime and iterations
times    = zeros(Nx, 1); % per state
maxIters = params.maxIters_;

% Cost matrix + constraint matrix
Cost   = build_mtx_cost(params);
Constr = build_mtx_constraint(sys, params);

nPhi = Nx*T + Nu*(T-1);

% ADMM variables
Phi = zeros(nPhi, Nx);

if nargin == 4 && ~isempty(warmStartIn) % warm start
    Psi    = warmStartIn.Psi_;
    Lambda = warmStartIn.Lambda_;
else
    Psi    = zeros(nPhi, Nx);
    Lambda = zeros(nPhi, Nx);
end

% Get sparsity of ADMM variables
PsiSupp  = get_sparsity_psi(sys, params);
PhiSupp  = PsiSupp(:, 1:Nx); % First block column
r   = assign_rows_phi(sys, T);
c   = assign_cols_phi(sys);
s_r = get_locality_row(PhiSupp);
s_c = get_locality_col(PhiSupp);

% Precalculate values for column update
[zabs, eyes, zabis] = precalculate_col_nominal(sys, T, s_c);

% Will be updated if adaptive ADMM is used
rho = params.rho_;

%% MPC
for iters=1:maxIters % ADMM loop
    Psi_prev = Psi;
    
    % Solve row-wise update for Phi
    Phi_rows = cell(nPhi, 1);
    for i = 1:Nx
        for row = r{i}
            x_loc = x0(s_r{row});
            cost  = Cost(row, row);

            if Constr(row, row) % has constraint
                if row <= T*Nx % state constraint
                    ub_ = params.stateUB_(i) / Constr(row, row);
                    lb_ = params.stateLB_(i) / Constr(row, row);
                else % input constraint
                    inputIdx = find(sys.B2(i,:));
                    ub_ = params.inputUB_(inputIdx) / Constr(row, row);
                    lb_ = params.inputLB_(inputIdx) / Constr(row, row);
                end
                ub = max(ub_,lb_); lb = min(ub_,lb_); % in case of negative signs
                solverMode = MPCSolverMode.Explicit;
            else % no constraint, use closed form
                solverMode = MPCSolverMode.ClosedForm;
            end
                        
            % User can override explicit solution and only use solver
            if params.useSolver_ && solverMode == MPCSolverMode.Explicit
                solverMode = MPCSolverMode.UseSolver;
            end
            
            if solverMode == MPCSolverMode.ClosedForm
                tic;
                Phi_rows{row} = row_nominal_closed(x_loc, Psi(row, s_r{row}), Lambda(row, s_r{row}), cost, rho);
                times(i) = times(i) + toc;
            elseif solverMode == MPCSolverMode.Explicit
                tic;
                Phi_rows{row} = row_nominal_explicit(x_loc, Psi(row, s_r{row}), Lambda(row, s_r{row}), ub, lb, cost, rho);
                times(i) = times(i) + toc;
            else
                [Phi_rows{row}, solverTime] = row_nominal_solver(x_loc, Psi(row, s_r{row}), Lambda(row, s_r{row}), ub, lb, cost, rho);
                times(i) = times(i) + solverTime;                    
            end
        end
    end
    
    % Build row-wise matrix Phi
    Phi = build_from_rows(r, s_r, Phi_rows, size(Phi));
        
    % Solve column-wise update for Psi
    Psi_cols = cell(Nx, 1);
    for i = 1:Nx
        tic;
        Psi_cols{i} = col_nominal_closed(Phi(s_c{i}, c{i}), Lambda(s_c{i}, c{i}), zabs{i}, eyes{i}, zabis{i});
        times(i) = times(i) + toc;        
    end
    
    % Build column-wise matrix Psi
    Psi = build_from_cols(c, s_c, Psi_cols, size(Psi));
    
    % Update Lambda matrix
    Lambda = Lambda + Phi - Psi;
    
    % Check convergence of ADMM
    converged = true;
    for i = 1:Nx
        phi_ = []; psi_ = []; psi_prev_ = [];
        for row = r{i}
            % Can't stack into small matrices; just use one big row
            phi_      = [phi_, Phi(row, s_r{row})];
            psi_      = [psi_, Psi(row, s_r{row})];
            psi_prev_ = [psi_prev_, Psi_prev(row, s_r{row})];
        end
        [conv, scale] = check_convergence(phi_, psi_, psi_, psi_prev_, params);

        if ~conv
            if params.has_adaptive_admm()
                rho = params.rho_;
            end
            converged = false; break; % if one fails, don't need to check the rest
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
u = Phi(1+Nx*T:Nx*T+Nu,:)*x0;
x = Phi(1+Nx:2*Nx,:)*x0; % since no noise, x_ref = x

% Running stats (runtime, iters)
stats            = MPCStats();
stats.time_      = mean(times);
stats.iters_     = iters;
stats.consIters_ = 0; % no consensus in this algorithm

% Warm start for next iteration
warmStartOut         = MPCWarmStart();
warmStartOut.Psi_    = Psi;
warmStartOut.Lambda_ = Lambda;

end