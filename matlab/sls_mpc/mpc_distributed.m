function [x, u, time, iters, consIters] = mpc_distributed(sys, x0, params)
% Inputs
%   sys     : LTISystem containing system matrices (A, B2) and Nx, Nu
%   x0      : Initial system state
%   params  : MPCParams containing parameters for mpc
% Outputs
%   x        : Next state if MPC input is used
%   u        : Current input as calculated by MPC
%   time     : Runtime (Steps 4+6+8+11) per state
%   iters    : ADMM iters per state (outer loop)
%   consIters: Average ADMM consensus iters per state, per iter (inner loop)
%              example: if outer loop ran 3 iters and inner loop ran
%                       2 iters, 5 iters, and 6 iters, then consIters = 4.3

%% Setup
sanity_check_actuation(sys)

if params.has_coupling()
    params.sanity_check_alg_2();
    sanity_check_coupling(sys, params);
else
    params.sanity_check_alg_1();
end

% For ease of notation
Nx = sys.Nx; Nu = sys.Nu; T = params.tFIR_;

maxIters     = params.maxIters_;
maxItersCons = params.maxItersCons_;

nPhi = Nx*T + Nu*(T-1);

% Track runtime and iterations
times        = zeros(Nx, 1); % per state
consIterList = zeros(maxIters,1);

% Get indices corresponding to rows / columns / localities
PsiSupp  = get_psi_sparsity(sys, params); % Toeplitz matrix
PhiSupp  = PsiSupp(:, 1:Nx);              % First block column
r   = assign_rows_phi(sys, T);
c   = assign_cols_phi(sys);
s_r = get_row_locality(PhiSupp);
s_c = get_col_locality(PhiSupp);

% Constraints and costs
C     = build_cost_mtx(params);
K     = build_constr_mtx(sys, params);

% ADMM variables
Phi    = zeros(nPhi, Nx);
Psi    = zeros(nPhi, Nx);
Lambda = zeros(nPhi, Nx);

% Coupling info and variables
cpIdx = get_coupling_indices_phi(C, K);
[rCp, rUcp, nValsCp] = sort_rows_coupled(r, cpIdx);
Ys = cell(nValsCp, 1);
Zs = cell(nValsCp, 1);

% Initialize Y and Z
for i=1:Nx
    for j = 1:length(rCp{i})
        row = rCp{i}{j};
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

    % Step 4: Solve for Phi
    Phi_rows = cell(nPhi, 1);
    
    % Solve for Phi for uncoupled rows
    for i = 1:Nx
        for row = rUcp{i}
            x_loc = x0(s_r{row}); % observe local state
            cost_ = C(row, row);

            if K(row, row) % has constraint
                if row <= T*Nx % state constraint
                    b1_ = params.stateUB_(i) / K(row, row);
                    b2_ = params.stateLB_(i) / K(row, row);
                else % input constraint
                    inputIdx = find(sys.B2(i,:));
                    b1_ = params.inputUB_(inputIdx) / K(row, row);
                    b2_ = params.inputLB_(inputIdx) / K(row, row);
                end
                b1  = max(b1_,b2_); b2 = min(b1_,b2_); % in case of negative signs
                solverMode = MPCSolverMode.Explicit;
            else % no constraint, use closed form
                solverMode = MPCSolverMode.ClosedForm;
            end
            
            % Terminal constraint, and row represents state at time T
            if params.terminalZeroConstr_ && row >= (T-1)*Nx + 1 && row <= T*Nx 
                b1 = 0; b2 = 0;
                solverMode = MPCSolverMode.Explicit;
            end
            
            if params.useSolver_ && solverMode == MPCSolverMode.Explicit
                solverMode = MPCSolverMode.UseSolver;
            end
            
            if solverMode == MPCSolverMode.ClosedForm
                tic;
                Phi_rows{row} = eqn_16a_closed(x_loc, Psi(row, s_r{row}), Lambda(row, s_r{row}), cost_, rho);
                times(i) = times(i) + toc;
            elseif solverMode == MPCSolverMode.Explicit
                tic;
                Phi_rows{row} = eqn_16a_explicit(x_loc, Psi(row, s_r{row}), Lambda(row, s_r{row}), b1, b2, cost_, rho);
                times(i) = times(i) + toc;
            else % use solver
                [Phi_rows{row}, solverTime] = eqn_16a_solver(x_loc, Psi(row, s_r{row}), Lambda(row, s_r{row}), b1, b2, cost_, rho);
                times(i) = times(i) + solverTime;                    
            end
        end
    end
    
    % Solve for Phi for coupled rows using ADMM consensus    
    if params.has_coupling()
        for consIter=1:maxItersCons % ADMM consensus (inner loop)
            Zs_prev = Zs;
            Xs      = cell(nValsCp, 1);

            for i = 1:Nx
                for row = rCp{i}
                    x_loc   = x0(s_r{row});     % observe local state
                    cps     = cpIdx{row};       % coupling indices for this row
                    selfIdx = find(cps == row); % index of "self-coupling" term

                    cost_ = C(row, cps);
                    k_    = K(row, cps); % constraint

                    if ~all(k_ == 0) % has constraint
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
                        [Phi_rows{row}, x_row] = eqn_20a_closed(x_loc, Psi(row, s_r{row}), Lambda(row, s_r{row}), ...
                                                 Ys{row}(cps), Zs(cps), cost_, selfIdx, params);
                        times(i) = times(i) + toc;
                    else % use solver
                        [Phi_rows{row}, x_row, solverTime] = eqn_20a_solver(x_loc, Psi(row, s_r{row}), Lambda(row, s_r{row}), ...
                                                             Ys{row}(cps), Zs(cps), cost_, k_, selfIdx, lb, ub, params);
                        times(i) = times(i) + solverTime;
                    end

                    Xs{row}             = zeros(nPhi, 1);
                    Xs{row}(cpIdx{row}) = x_row;
                end
            end

            % Step 6: Update Z (Step 5 implicitly done in this step)
            for i = 1:Nx
                tic;
                for row = rCp{i}
                    Zs{row} = 0;
                    for k = cpIdx{row}                                           
                        Zs{row} = Zs{row} + (Xs{k}(row)+Ys{k}{row})/length(cpIdx{row});
                    end
                end
                times(i) = times(i) + toc;
            end

            % Step 8: Update Y (Step 7 implicitly done in this step)
            for i = 1:Nx
                tic;
                for row = rCp{i}
                    for k = cpIdx{row}
                        Ys{row}{k} = Ys{row}{k} + Xs{row}(k) - Zs{k};
                    end
                end
                times(i) = times(i) + toc;                
            end

            % Step 9: Check convergence of ADMM consensus
            converged = true;
            for i = 1:Nx
                for row = rCp{i}
                    z_cp = zeros(nPhi, 1);
                    z_cp(cpIdx{row}) = [Zs{cpIdx{row}}];
                    
                    if ~check_convergence_cons(z_cp, Xs{row}, Zs{row}, Zs_prev{row}, params)
                        converged = false;
                        break; % if one fails, can stop checking the rest
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
    
    % Step 10: Build entire Phi matrix
    Phi = build_from_rows(r, s_r, Phi_rows, size(Phi));
        
    % Step 11: Solve (16b) to get local Psi
    Psi_cols = cell(Nx, 1);
    for i = 1:Nx
        tic;
        Psi_cols{i} = eqn_16b(Phi(s_c{i}, c{i}), Lambda(s_c{i}, c{i}), zabs{i}, eyes{i}, zabis{i});
        times(i) = times(i) + toc;        
    end
    
    % Step 12: Build entire Psi matrix
    Psi = build_from_cols(c, s_c, Psi_cols, size(Psi));
    
    % Step 13: Update Lambda
    Lambda = Lambda + Phi - Psi;
    
    % Step 14: Check convergence of ADMM (outer loop)
    converged = true;
    for i = 1:Nx
        phi_      = [];
        psi_      = [];
        psi_prev_ = [];
        for row = r{i}
            % Due to dimensionality issues, not stacking rows
            % Instead, just make one huge row
            % (since we're checking Frob norm, doesn't matter)
            phi_      = [phi_, Phi(row, s_r{row})];
            psi_      = [psi_, Psi(row, s_r{row})];
            psi_prev_ = [psi_prev_, Psi_prev(row, s_r{row})];
        end
        
        [conv, scale] = check_convergence(phi_, psi_, psi_, psi_prev_, params);
        
        if ~conv
            converged   = false;
            
            if params.has_adaptive_admm()
                % Alg2 functions use params instead of rho directly
                params.rho_ = min(params.rho_ * scale, params.rhoMax_);
                rho         = params.rho_;
            end
            
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
u = Phi(1+Nx*T:Nx*T+Nu,:)*x0;
x = Phi(1+Nx:2*Nx,:)*x0; % since no noise, x_ref = x

time      = mean(times);
consIters = mean(consIterList(1:iters));
params.rho_ = rho_original; % restore rho

end