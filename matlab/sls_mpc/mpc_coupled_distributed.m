function [x, u, stats, warmStartOut] = mpc_coupled_distributed(sys, x0, params, warmStartIn)
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
params.check_consensus_params();
params.sanity_check_dist();

% For ease of notation
Nx = sys.Nx; Nu = sys.Nu; T = params.tFIR_;

% Runtime and iterations
times        = zeros(Nx, 1); % per state
maxIters     = params.maxIters_;
maxItersCons = params.maxItersCons_;
consIterList = zeros(maxIters, 1); % for consensus

% Cost matrix
Cost   = build_mtx_cost(params);

% State/input/disturbance constraints
[H, h] = get_constraint_h(sys, params);

nH   = size(H, 1); % number of rows
nPhi = Nx*T + Nu*(T-1);

% ADMM variables
Phi   = zeros(nPhi, Nx);
Omega = zeros(nH, Nx);
    
if nargin == 4 && ~isempty(warmStartIn) % warm start
    fprintf('(Warm-started)\n');    
    Psi    = warmStartIn.Psi_;
    Lambda = warmStartIn.Lambda_;
else
    Psi    = zeros(nPhi, Nx);
    Lambda = zeros(nPhi + nH, Nx);
end

% Get sparsity of ADMM variables
HSupp = abs(H) > 0;
CSupp = abs(Cost) > 0;
PsiSupp = get_sparsity_psi(sys, params);
PsiSupp    = PsiSupp(:, 1:Nx); % first block column
PhiSupp    = CSupp * PsiSupp;

OmegaSupp  = zeros(0, Nx);
if nH > 0
    OmegaSupp  = abs(HSupp * PsiSupp) > 0; % since Omega = H * Phi
end
LambdaSupp = [PhiSupp; OmegaSupp];

% Which rows/cols of separable matrices are assigned to each subsystem
cPsi   = assign_cols_phi(sys);       % for Psi, Lambda
rPhi   = assign_rows_phi(sys, T);    % for Phi, Psi
rH_all = assign_rows_h(sys, params); % for Omega

% During row-wise update, can access neighbors' rows of Psi
% since we are not optimizing over Psi
rPsi_neighbors = get_local_neighbor_rows_phi(sys, params, rPhi);

s_rPsi    = get_locality_row(PsiSupp);
s_rPhi    = get_locality_row(PhiSupp);
s_rOmega  = get_locality_row(OmegaSupp);
s_rLambda = get_locality_row(LambdaSupp);

s_cPsi    = get_locality_col(PsiSupp);
s_cPhi    = get_locality_col(PhiSupp);
s_cOmega  = get_locality_col(OmegaSupp);
s_cLambda = get_locality_col(LambdaSupp);

% Precalculate values for column update
[cs, ms, gs] = precalculate_col_robust(sys, T, Cost, H, s_cPsi, s_cLambda);

% Will be updated if adaptive ADMM is used
rho = params.rho_;
mu  = params.mu_; % ease of notation

rH = rH_all;
% Terminal cost consensus
if params.terminal_cost_
    nHTerm    = length(params.terminal_h_);
    termHRows = nH-nHTerm+1:nH;
    
    % used in terminal consensus
    rHT       = cell(Nx, 1); % terminal rows, with consensus
    commsAdj  = sys.A ~= 0;
    neighbors = commsAdj^(params.locality_-1) ~= 0;
    
    ini = cell(Nx, 1); % neighbors, i.e. in_i   
    for i=1:Nx
        ini{i} = find(neighbors(i, :));        
        for rowH=rH{i}
            if ismember(rowH, termHRows)
                % move the row to rHT
                rHT{i}(end+1)      = rowH;
                rH{i}(rH{i}==rowH) = [];
            end
        end
    end
end
 
Y = zeros(Nx, 1);
Z = zeros(Nx, 1);
W = zeros(Nx, Nx);

%% MPC
for iters=1:maxIters % ADMM loop
    Psi_prev = Psi;

    % Solve row-wise update for Phi, Omega (no consensus)
    Phi_rows   = cell(nPhi, 1);
    Omega_rows = cell(nH, 1);    
    for i = 1:Nx     
        for row = rPhi{i}           
            [CPsiRow, s_rFull] = get_row_mp(Cost, Psi, s_rPsi, rPhi{i}, rPsi_neighbors{i}, row);
            x_loc = x0(s_rFull(s_rFull <= Nx)); 
            
            tic;
            Phi_rows{row} = row_nominal_closed(x_loc, CPsiRow, Lambda(row, s_rLambda{row}), 1, rho);
            times(i) = times(i) + toc;
        end
        
        for rowH = rH{i}
            rowLamb            = nPhi + rowH;
            [HPsiRow, s_rFull] = get_row_mp(H, Psi, s_rPsi, rPhi{i}, rPsi_neighbors{i}, rowH);
            x_loc = x0(s_rFull(s_rFull <= Nx)); 
            
            % User can override explicit solution and only use solver
            if params.useSolver_
                lenOmega = length(s_rOmega{rowH});
		        x_loc    = x0(s_rFull(s_rFull <= Nx)); 
		        a        = HPsiRow' - Lambda(rowLamb, s_rLambda{rowLamb})';
		                
		        C = zeros(lenOmega);
		        B = x_loc';
		        d = h(rowH);
		              
		        [Omega_rows{rowH}, solverTime] = row_lbn_solver(C, a, B, d);
                times(i) = times(i) + solverTime;
            else
                tic;
                Omega_rows{rowH} = row_nominal_explicit(x_loc, HPsiRow, Lambda(rowLamb, s_rLambda{rowLamb}), h(rowH), -inf, 0, rho);
                times(i) = times(i) + toc;
            end
        end
    end
    
    % Solve rows of Omega corresponding to terminal constraints via consensus
    if params.terminal_cost_
        for consIter=1:maxItersCons
            Z_prev = Z;
            X      = zeros(Nx, 1);
            
            for i=1:Nx
                C = []; B1 = []; B2 = []; a = [];
                
                for rowH=rHT{i}
                    rowLamb            = nPhi + rowH;
                    [HPsiRow, s_rFull] = get_row_mp(H, Psi, s_rPsi, rPhi{i}, rPsi_neighbors{i}, rowH);
                    x_loc              = x0(s_rFull(s_rFull <= Nx));

                    C  = blkdiag(C, rho/2 * eye(length(s_rOmega{rowH})));
                    B1 = blkdiag(B1, x_loc');
                    B2 = [B2; -h(rowH)];
                    a  = [a; rho/2 * (HPsiRow' - Lambda(rowLamb, s_rLambda{rowLamb})')];
                end
                C = blkdiag(C, mu/2);
                B = [B1                   B2;
                     zeros(2, size(B1,2)) [1; -1]];
                a = [a; mu/2 * Z(i)];
                F        = zeros(1, size(C,2));
                F(end)   = 1/Nx + Y(i);
                d        = zeros(size(B,1), 1);
                d(end-1) = 1;

                [v, solverTime] = row_consensus_solver(C, a, F, B, d);
                times(i)        = times(i) + solverTime;
                X(i)            = v(end); % consensus variable

                % Extract rows of Omega from answer
                rowStart = 1;
                for rowH=rHT{i}
                    rowEnd = rowStart + length(s_rOmega{rowH}) - 1;
                    Omega_rows{rowH} = v(rowStart:rowEnd);
                    rowStart = rowEnd + 1;
                end
            end
        
            % Update Z
            for i=1:Nx
                tic;
                loc   = ini{i}; % local patch/neighbors
                Z(i)  = sum(X(loc) + Z(loc) + (Y(loc) + W(i,loc)')/mu) / 2 / length(loc);
                times(i) = times(i) + toc;
            end
            
            % Update Y, W
            for i=1:Nx
                tic;
                loc      = ini{i};
                Y(i)     = Y(i) + mu*(X(i) - Z(i));
                W(i,loc) = W(i,loc) + mu*(Z(loc)' - Z(i)*ones(1, length(loc)));
                times(i) = times(i) + toc;
            end
        
            % Check convergence of ADMM consensus
            converged = true;
            for i = 1:Nx
                loc       = ini{i};
                primRes   = norm(X(i) - Z(i), 'fro') + norm(Z(loc) - Z(i)*ones(length(loc), 1), 'fro');
                dualRes   = norm(Z(i) - Z_prev(i), 'fro');
                converged = primRes <= params.eps_x_ && dualRes <= params.eps_z_;
                
                if ~converged
                    break; % exit the whole check loop
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
    
    % Build row-wise matrices Phi, Omega
    Phi   = build_from_rows(rPhi, s_rPhi, Phi_rows, size(Phi));    
    Omega = build_from_rows(rH_all, s_rOmega, Omega_rows, size(Omega));
    
    % Solve column-wise update for Psi
    Psi_cols = cell(Nx, 1);
    for i = 1:Nx
       v = [Phi(s_cPhi{i}, i); Omega(s_cOmega{i}, i)] + Lambda(s_cLambda{i}, i);
       tic;
       Psi_cols{i} = col_robust_closed(cs{i}, ms{i}, gs{i}, v);
       times(i) = times(i) + toc;
    end
    
    % Build column-wise matrix Psi
    Psi = build_from_cols(cPsi, s_cPsi, Psi_cols, size(Psi));    
    
    % Row-wise update for Lambda (split into two block rows)
    Lambda(1:nPhi, :)     = Lambda(1:nPhi, :) + Phi - Cost*Psi;
    if nH > 0
        Lambda(nPhi+1:end, :) = Lambda(nPhi+1:end, :) + Omega - H*Psi;        
    end
    
    % Check convergence of ADMM
    converged = true;
    for i = 1:Nx
        % Instead of stacking rows (dimension mismatch), make one huge row;
        % Doesn't matter for Frobenius norm anyway
        psi_      = [];
        psi_prev_ = [];
        prim1_    = []; % The I, H, Psi block
        prim2_    = []; % The Phi, Omega block
 
        for row = rPhi{i}        
            [CPsiRow, ~] = get_row_mp(Cost, Psi, s_rPsi, rPhi{i}, rPsi_neighbors{i}, row);
            psi_      = [psi_, Psi(row, s_rPsi{row})];
            psi_prev_ = [psi_prev_, Psi_prev(row, s_rPsi{row})];
             
            prim1_    = [prim1_, CPsiRow];
            prim2_    = [prim2_, Phi(row, s_rPhi{row})];                
        end
        for rowH = rH_all{i}
            [HPsiRow, ~] = get_row_mp(H, Psi, s_rPsi, rPhi{i}, rPsi_neighbors{i}, rowH);    
            prim1_ = [prim1_, HPsiRow];
            prim2_ = [prim2_, Omega(rowH, s_rOmega{rowH})];          
        end
        [conv, scale] = check_convergence(prim1_, prim2_, psi_, psi_prev_, params);
            
        if ~conv
            if params.has_adaptive_admm()
                rho = min(rho * scale, params.rhoMax_);
            end                
            converged = false; break; % if one fails, can stop checking the rest
        end
    end

    if converged
        break; % exit ADMM iterations
    end
end

if ~converged
    fprintf('ADMM reached %d iters and did not converge\n', maxIters);
end

% Compute control + next state
u = Psi(Nx*T+1:Nx*T+Nu, 1:Nx)*x0;
x = Psi(Nx+1:Nx*2, 1:Nx)*x0; % Next state, if no disturbance

% Running stats (runtime, iters)
stats            = MPCStats();
stats.time_      = mean(times); % average across all states
stats.iters_     = iters;
stats.consIters_ = mean(consIterList(1:iters));

% Warm start for next iteration
warmStartOut         = MPCWarmStart();
warmStartOut.Psi_    = Psi;
warmStartOut.Lambda_ = Lambda;

end