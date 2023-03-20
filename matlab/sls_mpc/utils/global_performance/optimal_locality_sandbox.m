clear all; clc;
warning off;

%% Grid example
seed          = 400;
gridSize      = 7;
T             = 15;
connectThresh = 0.6;
actDens       = 1.0;
Ts            = 0.2;

visualize = false;

params       = MPCParams();
params.tFIR_ = T+1; % Code and paper use different conventions

numNodes      = gridSize * gridSize; 
numActs       = round(actDens*numNodes);
rng(seed); % For actuated nodes
[adjMtx, nodeCoords, susceptMtx, inertiasInv, dampings] = generate_grid_topology(gridSize, connectThresh, seed);
actuatedNodes = randsample(numNodes, numActs);
sys = generate_grid_plant(actuatedNodes, adjMtx, susceptMtx, inertiasInv, dampings, Ts);

% Use custom communication structure for grid
sys.AComm = adjust_grid_sys_locality(sys.A);

%% Plot graph and actuated coordinates
if visualize
    plot_graph(adjMtx, nodeCoords, 'k')
    for i=1:numActs
        node = actuatedNodes(i);
        plot_special_vertex(node, nodeCoords, 'r')
    end
end

%% Get locality
tic;
locality = get_optimal_locality(sys, params);
toc

%% Setup for control problem
x0 = rand(sys.Nx, 1);

params.locality_ = locality;
params.QSqrt_    = diag(rand(sys.Nx, 1));
params.RSqrt_    = diag(rand(sys.Nu, 1));

Nx = sys.Nx; Nu = sys.Nu; T = params.tFIR_;
nPhi  = Nx*T + Nu*(T-1);
IO    = [eye(Nx); zeros(Nx*(T-1), Nx)];
QSqrt = params.QSqrt_;
RSqrt = params.RSqrt_;
ZAB   = get_constraint_zab(sys, T);

%% Calculate global trajectory
fprintf('Calculating global trajectory\n');

cvx_begin quiet
variable xs(Nx, T)
variable us(Nu, T-1)

xs(:,1) == x0; % Initial condition
obj1    = 0; % Objective function
for t = 1:T-1
    obj1 = obj1 + xs(:,t)'*QSqrt*QSqrt*xs(:,t) ...
                + us(:,t)'*RSqrt*RSqrt*us(:,t);
            
    % Dynamics
    xs(:,t+1) == sys.A*xs(:,t) + sys.B2*us(:,t);
            
end
obj1 = obj1 + + xs(:,T)'*QSqrt*QSqrt*xs(:,T); % Terminal cost

minimize(obj1)
cvx_end

%% Calculate localized trajectory [Adapted from mpc_centralized]
fprintf('Calculating localized trajectory\n');

PsiSupp = get_sparsity_psi(sys, params); 
PsiSupp = PsiSupp(:, 1:Nx);
suppSizePsi = sum(sum(PsiSupp));

cvx_begin quiet
expression Psi(nPhi, Nx)
variable PsiSuppVals(suppSizePsi, 1)
Psi(PsiSupp) = PsiSuppVals; 

% Set up objective function
obj2 = 0;
for k = 1:T
    kx         = get_range(k, Nx); % rows of Psi representing Psix
    vect       = QSqrt*Psi(kx, 1:Nx)*x0;
    obj2 = obj2 + vect'*vect;
end
for k = 1:T-1
    ku         = Nx*T + get_range(k, Nu); % rows of Psi representing Psiu
    vect       = RSqrt*Psi(ku, 1:Nx)*x0;
    obj2 = obj2 + vect'*vect;
end

% Dynamics constraints
% Note: have to use this formulation b/c otherwise  cvx reports infeasible
%       for ZAB*Psi == IO, even if we set low precision
EPS = 1e-8;
norm(ZAB*Psi - IO, 'fro') <= EPS;
minimize(obj2)
cvx_end

trajDiff = norm([vec(xs); vec(us)] - Psi*x0) / norm(xs);
objDiff  = norm(obj2 - obj1) / obj1;

fprintf('Diff between global and localized objectives: %.3e\n', objDiff);
fprintf('Diff between global and localized trajectory: %.3e\n', trajDiff);


