function [x, u, time, iters, consIters] = mpc_distributed(sys, x0, params)
% Inputs
%   sys     : LTISystem containing system matrices (A, B2) and Nx, Nu
%   x0      : Initial system state
%   params  : MPCParams containing parameters for mpc
% Outputs
%   x        : Next state if MPC input is used
%   u        : Current input as calculated by MPC
%   time     : Runtime per subsystem; time taken to calculate
%              Phi (row update), Psi (col update), and consensus
%              updates on X and Z if applicable
%   iters    : ADMM iters per state (outer loop)
%   consIters: Average ADMM consensus iters per state, per iter (inner loop)
%              example: if outer loop ran 3 iters and inner loop ran
%                       2 iters, 5 iters, and 6 iters, then consIters = 4.3

%% Setup
sanity_check_actuation(sys);
sanity_check_coupling(sys, params);
params.sanity_check_dist();

% For ease of notation
Nx = sys.Nx; Nu = sys.Nu; T = params.tFIR_;

maxIters     = params.maxIters_;
maxItersCons = params.maxItersCons_;

nPhi = Nx*T + Nu*(T-1);

% Track runtime and iterations
times        = zeros(Nx, 1); % per state
consIterList = zeros(maxIters,1);

% Get indices corresponding to rows / columns / localities
PsiSupp  = get_psi_sparsity(sys, params);
PhiSupp  = PsiSupp(:, 1:Nx); % First block column
r   = assign_rows_phi(sys, T);
c   = assign_cols_phi(sys);
s_r = get_row_locality(PhiSupp);
s_c = get_col_locality(PhiSupp);

% Constraints and costs
Cost   = build_cost_mtx(params);
Constr = build_constr_mtx(sys, params);

% ADMM variables
Phi    = zeros(nPhi, Nx);
Psi    = zeros(nPhi, Nx);
Lambda = zeros(nPhi, Nx);

% Coupling info and variables
cpIdx = get_coupling_indices_phi(Cost, Constr);
[rCp, rUcp, nValsCp] = sort_rows_coupled(r, cpIdx);
Ys = cell(nValsCp, 1);
Zs = cell(nValsCp, 1);

% Initialize Y and Z
for i=1:Nx
    for row = rCp{i}
        Zs{row} = 0;
        for k = cpIdx{row}
            Ys{row}{k} = 0;
        end
    end
end

% Precalculate items for column-wise update
[zabs, eyes, zabis] = precalculate_col(sys, T, s_c);

% We will change params.rho_ but we want to set it back afterwards
rho_original = params.rho_;
rho          = params.rho_;

%% MPC
for iters=1:maxIters % ADMM (outer loop)
    Psi_prev = Psi;
    Phi_rows = cell(nPhi, 1);
    
    % Solve for uncoupled rows in Phi
    for i = 1:Nx
        for row = rUcp{i}
            x_loc = x0(s_r{row}); % observe local state
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
            
            % Terminal constraint, and row represents state at time T
            if params.terminalZeroConstr_ && row >= (T-1)*Nx + 1 && row <= T*Nx 
                ub = 0; lb = 0;
                solverMode = MPCSolverMode.Explicit;
            end
            
            % User can choose to override explicit mode and only use solver
            if params.useSolver_ && solverMode == MPCSolverMode.Explicit
                solverMode = MPCSolverMode.UseSolver;
            end
            
            if solverMode == MPCSolverMode.ClosedForm
                tic;
                Phi_rows{row} = mpc_row_closed(x_loc, Psi(row, s_r{row}), Lambda(row, s_r{row}), cost, rho);
                times(i) = times(i) + toc;
            elseif solverMode == MPCSolverMode.Explicit
                tic;
                Phi_rows{row} = mpc_row_explicit(x_loc, Psi(row, s_r{row}), Lambda(row, s_r{row}), ub, lb, cost, rho);
                times(i) = times(i) + toc;
            else % use solver
                [Phi_rows{row}, solverTime] = mpc_row_solver(x_loc, Psi(row, s_r{row}), Lambda(row, s_r{row}), ub, lb, cost, rho);
                times(i) = times(i) + solverTime;                    
            end
        end
    end
    
    % Solve for coupled rows in Phi using ADMM consensus    
    if params.has_coupling()
        for consIter=1:maxItersCons % ADMM consensus (inner loop)
            Zs_prev = Zs;
            Xs      = cell(nValsCp, 1);

            for i = 1:Nx
                for row = rCp{i}
                    x_loc = x0(s_r{row});     % observe local state
                    cp    = cpIdx{row};       % coupling indices for this row
                    sIdx  = find(cp == row);  % index of "self-coupling" term
                    cost   = Cost(row, cp);
                    constr = Constr(row, cp);
                    
                    if ~all(constr == 0) % has constraint
                        if row <= Nx*T % is state
                            lb  = params.stateLB_(i);
                            ub  = params.stateUB_(i);
                        else % is input
                            inputIdx = find(sys.B2(i,:));
                            lb = params.inputLB_(inputIdx);
                            ub = params.inputUB_(inputIdx);
                        end
                        solverMode = MPCSolverMode.UseSolver;
                    else % no constraint, use closed form
                        solverMode = MPCSolverMode.ClosedForm;
                    end

                    if solverMode == MPCSolverMode.ClosedForm
                        tic;
                        [Phi_rows{row}, x_] = mpc_coupled_row_closed(x_loc, Psi(row, s_r{row}), Lambda(row, s_r{row}), ...
                                                 Ys{row}(cp), Zs(cp), cost, sIdx, params);
                        times(i) = times(i) + toc;
                    else % use solver
                        [Phi_rows{row}, x_, solverTime] = mpc_coupled_row_solver(x_loc, Psi(row, s_r{row}), Lambda(row, s_r{row}), ...
                                                             Ys{row}(cp), Zs(cp), cost, constr, sIdx, lb, ub, params);
                        times(i) = times(i) + solverTime;
                    end

                    Xs{row}             = zeros(nPhi, 1);
                    Xs{row}(cpIdx{row}) = x_;
                end
            end

            % Update Z for consensus
            for i = 1:Nx
                tic; 
                Zs = cons_z_update(Xs, Ys, Zs, rCp{i}, cpIdx);
                times(i) = times(i) + toc;
            end

            % Update Y for consensus
            for i = 1:Nx
                Ys = cons_y_update(Xs, Ys, Zs, rCp{i}, cpIdx);
            end

            % Check convergence of ADMM consensus
            converged = true;
            for i = 1:Nx
                for row = rCp{i}
                    z_cp = zeros(nPhi, 1);
                    z_cp(cpIdx{row}) = [Zs{cpIdx{row}}];
                    
                    if ~check_convergence_cons(z_cp, Xs{row}, Zs{row}, Zs_prev{row}, params)
                        converged = false; break; % if one fails, don't need to check the rest
                    end
                end
            end

            if converged
                break; % exit ADMM consensus iterations
            end
        end
        
        if ~converged
            fprintf('ADMM consensus reached %d iters and did not converge\n', maxItersCons);
        end
        
        consIterList(iters) = consIter;
    end 
    
    % Build Phi matrix
    Phi = build_from_rows(r, s_r, Phi_rows, size(Phi));
        
    % Solve for columns of Psi
    Psi_cols = cell(Nx, 1);
    for i = 1:Nx
        tic;
        Psi_cols{i} = mpc_col_closed(Phi(s_c{i}, c{i}), Lambda(s_c{i}, c{i}), zabs{i}, eyes{i}, zabis{i});
        times(i) = times(i) + toc;        
    end
    
    % Build Psi matrix
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
                % Coupled functions use params instead of rho directly
                params.rho_ = min(params.rho_ * scale, params.rhoMax_);
                rho         = params.rho_;
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

time      = mean(times);
consIters = mean(consIterList(1:iters));
params.rho_ = rho_original; % restore rho

end