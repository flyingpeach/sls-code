function [x, u, stats, warmStartOut] = rmpc_distributed(sys, x0, params, warmStartIn)
% Inputs
%   sys     : LTISystem containing system matrices (A, B2) and Nx, Nu 
%   x0      : Starting system state
%   params  : MPCParams containing parameters for mpc
%   (optional): warmStartIn: MPCWarmStart containing Psi, Lambda from
%               previous timestep
% Outputs
%   x       : Next state if MPC input is used and no disturbance
%   u       : Current input as calculated by MPC
%   stats    : MPCStats containing runtime/iteration stats
%              time: Runtime per subsystem; (row + column update)
%              iters: ADMM iters per state (outer loop)
%   warmStartOut: MPCWarmStart containing Psi, Lambda to be used for
%                 next time iteration, if applicable

%% Setup
sanity_check_actuation(sys);
sanity_check_coupling(sys, params);
params.sanity_check_dist();

if ~params.has_state_cons() && ~params.has_input_cons()
    mpc_error('No state/input constraints were specified for Robust MPC. This is not supported. \nIf you are using sls_mpc.m, clear the disturbance params and try again.');
end

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
Phi = zeros(nPhi, Nx);
if nargin == 4 && ~isempty(warmStartIn) % warm start
    fprintf('(Warm-started)\n');
    Psi    = warmStartIn.Psi_;
    Lambda = warmStartIn.Lambda_;
else
    Psi    = zeros(nPhi, Nx*T);
    Lambda = zeros(nPhi + nH, Nx*T);
end

% Get sparsity of ADMM variables
zBlock   = zeros(nPhi, Nx*(T-1));
HSupp    = abs(H) > 0;
CSupp    = abs(Cost) > 0;
PsiSupp  = get_sparsity_psi(sys, params);
Psi1Supp = PsiSupp(:, 1:Nx); % first block column
PhiSupp  = CSupp * Psi1Supp;

% Which rows/cols of separable matrices are assigned to each subsystem
cPsi = assign_cols_psi(sys, T);      % for Psi, Lambda
rPhi = assign_rows_phi(sys, T);      % for Phi, Psi
rH_all = assign_rows_h(sys, params); % for Omega, Xi

% During row-wise update, can access neighbors' rows of Psi
% since we are not optimizing over Psi
rPsi_neighbors = get_local_neighbor_rows_phi(sys, params, rPhi);

if params.has_polytopic_noise()
    [G, g] = get_constraint_g(params);
    nG     = size(G, 1);
    Omega  = zeros(nH, Nx);
    Xi     = zeros(nH, nG);

    GSupp      = abs(G) > 0;
    OmegaSupp  = abs(HSupp * Psi1Supp) > 0; % since Omega = H * Phi
    XiSupp     = get_sparsity_xi(H, G, Nx, PsiSupp);
    XiGSupp    = abs(XiSupp * GSupp) > 0; % sparsity of Xi * G
    LambdaSupp = [PhiSupp zBlock; OmegaSupp XiGSupp];
    
    cXi       = assign_cols_xi(sys, T); 
    s_rXi     = get_locality_row(XiSupp);
    s_cXi     = get_locality_col(XiSupp);
    % Pre-process G for easy multiplication with local Xi
    gts = precalculate_g_trimmed(G, s_rXi, nH);

else % locally bounded noise
    sigma = params.locNoiseBound_;
    Omega = zeros(nH, Nx*T);
    
    OmegaSupp  = abs(HSupp * PsiSupp) > 0; % Omega = H * Psi
    LambdaSupp = [PhiSupp zBlock; OmegaSupp]; 
end

s_rPsi    = get_locality_row(PsiSupp);
s_rPsi1   = get_locality_row(Psi1Supp);
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

rH = rH_all;
% Terminal cost consensus
if params.terminal_cost_
    mu        = params.mu_; % ease of notation
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
    if ~mod(iters, 250)
        fprintf('rMPC Iteration: %d\n', iters);
    end
    Psi_prev = Psi;

    % Solve row-wise update for Phi, Omega, Xi (no consensus)
    Phi_rows   = cell(nPhi, 1);
    Omega_rows = cell(nH, 1);
    Xi_rows    = cell(nH, 1);
    for i = 1:Nx     
        for row = rPhi{i}           
            [CPsiRow, s_rFull] = get_row_mp(Cost, Psi, s_rPsi1, rPhi{i}, rPsi_neighbors{i}, row);
            x_loc = x0(s_rFull(s_rFull <= Nx)); 
            
            tic;
            Phi_rows{row} = row_nominal_closed(x_loc, CPsiRow, Lambda(row, s_rLambda{row}), 1, rho);
            times(i) = times(i) + toc;
        end
        
        for rowH = rH{i}
            rowLamb            = nPhi + rowH;
            [HPsiRow, s_rFull] = get_row_mp(H, Psi, s_rPsi, rPhi{i}, rPsi_neighbors{i}, rowH);
            lenOmega = length(s_rOmega{rowH});
            x_loc    = x0(s_rFull(s_rFull <= Nx)); 
            a        = HPsiRow' - Lambda(rowLamb, s_rLambda{rowLamb})';
                
            if params.has_polytopic_noise()
                lenXi    = length(s_rXi{rowH});
                if lenXi == 0 % Xi is empty here, so it's just Omega
                    C = eye(lenOmega);
                    B = x_loc';
                    d = h(rowH);
                else                    
                    C  = blkdiag(eye(lenOmega), gts{rowH}');                     
                    g_ = g(s_rXi{rowH});                  
                    B  = [x_loc' g_'; zeros(lenXi, lenOmega) -eye(lenXi)];
                    d  = [h(rowH); zeros(lenXi, 1)];
                end
                F = zeros(1, size(C,2));
                
                [omegaxi, solverTime] = row_general_solver(C, a, F, B, d);
                Omega_rows{rowH} = omegaxi(1:lenOmega);                
                if lenXi ~= 0
                    Xi_rows{rowH} = omegaxi(lenOmega+1:end);
                end
            else % locally bounded noise
                lenOmega2 = length(s_rOmega{rowH}(s_rOmega{rowH} > Nx));
                lenOmega1 = lenOmega - lenOmega2;
                
                C = sigma * blkdiag(zeros(lenOmega1), eye(lenOmega2));
                B = [x_loc' zeros(1, lenOmega2)];
                d = h(rowH);
                
                [Omega_rows{rowH}, solverTime] = row_lbn_solver(C, a, B, d);
            end            
            times(i) = times(i) + solverTime;
         end
    end
    
    % Solve rows of Omega, Xi corresponding to terminal constraints via consensus
    if params.terminal_cost_
        for consIter=1:maxItersCons
            Z_prev = Z;
            X      = zeros(Nx, 1);
            
            for i=1:Nx
                C = []; B1 = []; B2 = []; B3 = []; a = [];
                
                for rowH=rHT{i}
                    lenOm = length(s_rOmega{rowH}); lenXi = length(s_rXi{rowH});
                    rowLamb            = nPhi + rowH;
                    [HPsiRow, s_rFull] = get_row_mp(H, Psi, s_rPsi, rPhi{i}, rPsi_neighbors{i}, rowH);
                    x_loc              = x0(s_rFull(s_rFull <= Nx));

                    C  = blkdiag(C, rho/2 * eye(lenOm), rho/2 * gts{rowH}');
                    B1 = blkdiag(B1, [x_loc' g(s_rXi{rowH})']); 
                    B2 = [B2; -h(rowH)];
                    B3 = blkdiag(B3, [zeros(lenXi, lenOm) -eye(lenXi)]);
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
                
                B = [B; B3 zeros(size(B3,1), 1)]; 
                d = [d; zeros(size(B3,1), 1)]; % non-negative Xi

                [v, solverTime] = row_general_solver(C, a, F, B, d);
                times(i)        = times(i) + solverTime;
                X(i)            = v(end); % consensus variable

                % Extract rows of Omega, Xi from answer
                rowStart = 1;
                for rowH=rHT{i}
                    rowEnd = rowStart + length(s_rOmega{rowH}) - 1;
                    Omega_rows{rowH} = v(rowStart:rowEnd);
                    rowStart = rowEnd + 1;
                    rowEnd = rowStart + length(s_rXi{rowH}) - 1;
                    Xi_rows{rowH} = v(rowStart:rowEnd);
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
        
    if params.has_polytopic_noise()
        Xi = build_from_rows(rH_all, s_rXi, Xi_rows, size(Xi));    
    end
    
    % Solve column-wise update for Psi
    Psi_cols = cell(Nx*T, 1);
    for i = 1:Nx
        for col = cPsi{i}
            if col <= Nx
                v = [Phi(s_cPhi{col}, col); Omega(s_cOmega{col}, col)] + Lambda(s_cLambda{col}, col);
            else
                if params.has_polytopic_noise()
                    colG   = col - Nx;
                    XiGCol = get_col_xig(G, Xi, s_cXi, cXi{i}, colG);                
                    v      = XiGCol + Lambda(s_cLambda{col}, col);
                else % locally bounded noise
                    v = Omega(s_cOmega{col}, col) + Lambda(s_cLambda{col}, col);
                end
            end
            tic;
            Psi_cols{col} = col_robust_closed(cs{col}, ms{col}, gs{col}, v);
            times(i) = times(i) + toc;
        end
    end

    % Build column-wise matrix Psi
    Psi = build_from_cols(cPsi, s_cPsi, Psi_cols, size(Psi));    
                     
    % Row-wise update for Lambda (split into two block rows)
    Lambda(1:nPhi, :)     = Lambda(1:nPhi, :) + [Phi zBlock] - [Cost*Psi(:, 1:Nx) zBlock];
    
    if params.has_polytopic_noise()
        Lambda(nPhi+1:end, :) = Lambda(nPhi+1:end, :) + [Omega Xi*G] - H*Psi;
    else % locally bounded noise
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
            [CPsiRow, ~] = get_row_mp(Cost, Psi, s_rPsi1, rPhi{i}, rPsi_neighbors{i}, row);
            psi_      = [psi_, Psi(row, s_rPsi{row})];
            psi_prev_ = [psi_prev_, Psi_prev(row, s_rPsi{row})];
             
            prim1_    = [prim1_, CPsiRow];
            prim2_    = [prim2_, Phi(row, s_rPhi{row})];                
        end
        for rowH = rH_all{i}
            [HPsiRow, ~] = get_row_mp(H, Psi, s_rPsi, rPhi{i}, rPsi_neighbors{i}, rowH);    
            prim1_    = [prim1_, HPsiRow];

            prim2_ = [prim2_, Omega(rowH, s_rOmega{rowH})];

            if params.has_polytopic_noise()
                XiGRow = Xi(rowH, s_rXi{rowH}) * gts{rowH};
                prim2_ = [prim2_, XiGRow];
            end                
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