function [x, u, time] = mpc_centralized(sys, x0, params)
% Note: this is a per-timestep algorithm
params.sanity_check_cent();

Nx = sys.Nx; Nu = sys.Nu; T  = params.tFIR_;   

QSqrt = params.QSqrt_;
RSqrt = params.RSqrt_; % note this is not related to SLS R

% Record time of inner solver instead of full cvx solver
% Note: this code only works if we're using gurobi as a solver
% For debugging: find cvx_run_solver and check varargout
FILENAME = 'gData.mat';
cvx_solver_settings('dumpfile', FILENAME);

% Constraint matrices
ZAB     = get_constraint_zab(sys, T);
[H, h]  = get_constraint_h(sys, params);
PsiSupp = get_sparsity_psi(sys, params);
nH      = size(H, 1);

suppSizePsi = sum(sum(PsiSupp));

if params.has_polytopic_noise()
    [G, g] = get_constraint_g(params);
    XiSupp  = get_sparsity_xi(H, G, Nx, PsiSupp);
    suppSizeXi  = sum(sum(XiSupp));
end

cvx_begin quiet

expression Psi(Nx*T + Nu*(T-1), Nx*T)
variable PsiSuppVals(suppSizePsi)
Psi(PsiSupp) = PsiSuppVals;

if params.has_polytopic_noise()
    expression Xi(nH, size(G,1))
    variable XiSuppVals(suppSizeXi) nonnegative
    Xi(XiSupp)   = XiSuppVals;
end

% Set up objective function
objective = 0;
for k = 1:T
    kx        = get_range(k, Nx); % rows of Psi representing Psix
    state     = Psi(kx, 1:Nx)*x0;
    vect      = vec(QSqrt*state);
    objective = objective + vect'*vect;
    end
for k = 1:T-1
    ku        = Nx*T + get_range(k, Nu); % rows of Psi representing Psiu
    input     = Psi(ku, 1:Nx)*x0;
    vect      = vec(RSqrt*input);
    objective = objective + vect'*vect;
end

% Constraints
ZAB*Psi == eye(Nx*T); % dynamics constraints

if params.is_robust()
    if params.has_polytopic_noise()
        H*Psi(:, 1:Nx)*x0 + Xi*g <= h;
        H*Psi(:, Nx+1:end)       == Xi*G;
    else % locally bounded noise
        sigma = params.locNoiseBound_;
        for rowH = 1:nH
            H(rowH, :)*Psi(:, 1:Nx)*x0 + sigma*norm(H(rowH, :)*Psi(:, Nx+1:end)) <= h(rowH);
        end
    end
else % no noise expected
    if nH > 0
        H*Psi(:, 1:Nx)*x0 <= h;
    end
end

% don't need to enforce <= 1 because it's already done by constraints
if params.terminal_cost_
    variable eta nonnegative 
    if params.has_polytopic_noise()
        nHTerm  = length(params.terminal_h_);
        termIdx = nH-nHTerm+1:nH;       
        params.terminal_H_*Psi(Nx*(T-1)+1:Nx*T, 1:Nx)*x0 + Xi(termIdx,:)*g <= eta*params.terminal_h_;        
    else % nominal (since local noise not supported)
        params.terminal_H_*Psi(Nx*(T-1)+1:Nx*T, 1:Nx)*x0 <= eta*params.terminal_h_;
    end
    objective = objective + eta;
end

minimize(objective)
cvx_end

% Note: this code only works if we're using gurobi as a solver
gData = load(FILENAME);
time  = gData.res.runtime;

% Compute control + state
u = Psi(Nx*T+1:Nx*T+Nu, 1:Nx)*x0; % First block matrix of Psiu * x0
x = Psi(Nx+1:2*Nx, 1:Nx)*x0;      % Second block matrix of Psix * x0
                                  % Assumes no disturbance
end