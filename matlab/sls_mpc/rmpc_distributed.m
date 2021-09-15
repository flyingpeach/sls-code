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
times    = zeros(Nx, 1); % per state
maxIters = params.maxIters_;

% Cost matrix
Cost   = build_cost_mtx(params);

% State/input/disturbance constraints
[H, h] = get_sys_constraints(sys, params);

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
zBlock = zeros(nPhi, Nx*(T-1));
HSupp = abs(H) > 0;
CSupp = abs(Cost) > 0;
PsiSupp   = get_psi_sparsity(sys, params);
Psi1Supp  = PsiSupp(:, 1:Nx); % first block column
PhiSupp   = CSupp * Psi1Supp;

% Which rows/cols of separable matrices are assigned to each subsystem
cPsi = assign_cols_psi(sys, T);      % for Psi, Lambda
rPhi = assign_rows_phi(sys, T);      % for Phi, Psi
rH   = assign_rows_h(sys, params);   % for Omega, Xi

% During row-wise update, can access neighbors' rows of Psi
% since we are not optimizing over Psi
rPsi_neighbors = get_local_neighbor_rows_phi(sys, params, rPhi);

if params.has_polytopic_noise()
    [G, g] = get_dist_constraints(params);
    nG     = size(G, 1);
    Omega  = zeros(nH, Nx);
    Xi     = zeros(nH, nG);

    GSupp      = abs(G) > 0;
    OmegaSupp  = abs(HSupp * Psi1Supp) > 0; % since Omega = H * Phi
    XiSupp     = get_xi_sparsity(PsiSupp, H, G, Nx);
    XiGSupp    = abs(XiSupp * GSupp) > 0; % sparsity of Xi * G
    LambdaSupp = [PhiSupp zBlock; OmegaSupp XiGSupp];
    
    cXi       = assign_cols_xi(sys, T); 
    s_rXi     = get_row_locality(XiSupp);
    s_cXi     = get_col_locality(XiSupp);
    % Pre-process G for easy multiplication with local Xi
    gts = precalculate_g_trimmed(G, s_rXi, nH);

else % locally bounded noise
    sigma = params.locNoiseBound_;
    Omega = zeros(nH, Nx*T);
    
    OmegaSupp  = abs(HSupp * PsiSupp) > 0; % Omega = H * Psi
    LambdaSupp = [PhiSupp zBlock; OmegaSupp]; 
end

s_rPsi    = get_row_locality(PsiSupp);
s_rPsi1   = get_row_locality(Psi1Supp);
s_rPhi    = get_row_locality(PhiSupp);
s_rOmega  = get_row_locality(OmegaSupp);
s_rLambda = get_row_locality(LambdaSupp);

s_cPsi    = get_col_locality(PsiSupp);
s_cPhi    = get_col_locality(PhiSupp);
s_cOmega  = get_col_locality(OmegaSupp);
s_cLambda = get_col_locality(LambdaSupp);

% Precalculate values for column update
[cs, ms, gs] = precalculate_robust_col(sys, T, Cost, H, s_cPsi, s_cLambda);

% Will be updated if adaptive ADMM is used
rho = params.rho_;

%% MPC
for iters=1:maxIters % ADMM loop
    if ~mod(iters, 250)
        fprintf('rMPC Iteration: %d\n', iters);
    end
    Psi_prev = Psi;

    % Solve row-wise update for Phi, Omega, Xi
    Phi_rows   = cell(nPhi, 1);
    Omega_rows = cell(nH, 1);
    Xi_rows    = cell(nH, 1);
    for i = 1:Nx     
        for row = rPhi{i}           
            [CPsiRow, s_rFull] = get_cpsi_row(Cost, Psi, s_rPsi1, rPhi{i}, rPsi_neighbors{i}, row);
            x_loc = x0(s_rFull(s_rFull <= Nx)); 
            
            tic;
            Phi_rows{row} = mpc_row_closed(x_loc, CPsiRow, Lambda(row, s_rLambda{row}), rho);
            times(i) = times(i) + toc;
        end
        
        for rowH = rH{i}
            rowLamb            = nPhi + rowH;
            [HPsiRow, s_rFull] = get_hp_row(H, Psi, s_rPsi, rPhi{i}, rPsi_neighbors{i}, rowH);
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

                [omegaxi, solverTime] = calc_robust_row_solver(C, a, B, d);
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
                
                [Omega_rows{rowH}, solverTime] = calc_lbn_row_solver(C, a, B, d);
            end            
            times(i) = times(i) + solverTime;
         end
    end
    
    % Build row-wise matrices Phi, Omega, Xi
    Phi   = build_from_rows(rPhi, s_rPhi, Phi_rows, size(Phi));
    Omega = build_from_rows(rH, s_rOmega, Omega_rows, size(Omega));
    
    if params.has_polytopic_noise()
        Xi = build_from_rows(rH, s_rXi, Xi_rows, size(Xi));    
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
                    XiGCol = get_xig_col(G, Xi, s_cXi, cXi{i}, colG);                
                    v      = XiGCol + Lambda(s_cLambda{col}, col);
                else % locally bounded noise
                    v = Omega(s_cOmega{col}, col) + Lambda(s_cLambda{col}, col);
                end
            end
            tic;
            Psi_cols{col} = calc_robust_col_closed(cs{col}, ms{col}, gs{col}, v);
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
            [CPsiRow, ~] = get_cpsi_row(Cost, Psi, s_rPsi1, rPhi{i}, rPsi_neighbors{i}, row);
            psi_      = [psi_, Psi(row, s_rPsi{row})];
            psi_prev_ = [psi_prev_, Psi_prev(row, s_rPsi{row})];
             
            prim1_    = [prim1_, CPsiRow];
            prim2_    = [prim2_, Phi(row, s_rPhi{row})];                
        end
        for rowH = rH{i}
            [HPsiRow, ~] = get_hp_row(H, Psi, s_rPsi, rPhi{i}, rPsi_neighbors{i}, rowH);    
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
stats.consIters_ = 0; % no consensus in this algorithm

% Warm start for next iteration
warmStartOut         = MPCWarmStart();
warmStartOut.Psi_    = Psi;
warmStartOut.Lambda_ = Lambda;

end