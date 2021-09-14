function [x, u, time] = rmpc_centralized(sys, x0, params)
% Note: this is a per-timestep algorithm
params.sanity_check_cent();

Nx = sys.Nx; Nu = sys.Nu; T  = params.tFIR_;   

QSqrt = params.QSqrt_;
RSqrt = params.RSqrt_; % note this is not related to SLS R

time = 0;

% Constraint matrices
ZAB     = get_sls_constraint(sys, T);
[H, h]  = get_sys_constraints(sys, params);
PsiSupp = get_psi_sparsity(sys, params);
nH      = size(H, 1);

suppSizePsi = sum(sum(PsiSupp));

if params.has_polytopic_noise()
    [G, g] = get_dist_constraints(params);
    XiSupp  = get_xi_sparsity(PsiSupp, H, G, Nx);
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

tic;
minimize(objective)
subject to
ZAB*Psi == eye(Nx*T); % dynamics constraints

if params.has_polytopic_noise()
    H*Psi(:, 1:Nx)*x0 + Xi*g <= h;
    H*Psi(:, Nx+1:end)       == Xi*G;
else % locally bounded noise
    sigma = params.locNoiseBound_;
    for rowH = 1:nH
        H(rowH, :)*Psi(:, 1:Nx)*x0 + sigma*norm(H(rowH, :)*Psi(:, Nx+1:end)) <= h(rowH);
    end
end

time = time + toc;
cvx_end

% Compute control + state
u = Psi(Nx*T+1:Nx*T+Nu, 1:Nx)*x0; % First block matrix of Psiu * x0
x = Psi(Nx+1:2*Nx, 1:Nx)*x0;      % Second block matrix of Psix * x0
                                  % Assumes no disturbance
end