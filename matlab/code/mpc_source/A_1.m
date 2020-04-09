% Algorithm I, System A

%% Setup
setup_system_a;
setup_sls_constr;
setup_loc_constr;

% Weights
Q = eye(Nx);
S = diag(ones(Nu,1));

%% Distributed MPC
rho   = 5; % admm parameter
eps_d = 1e-3; % convergence criterion for ||Psi(k+1) - Psi(k)||
eps_p = 1e-4; % convergence criterion for ||Phi(k+1) - Psi(k+1)||
maxIters = 5000;

[x, u] = alg1_fn(Nx, Nu, A, B, d, tFIR, tSim, x0, ...
                 eps_d, eps_p, rho, maxIters);

%% Centralized MPC (for validation + comparison)
[xVal, uVal] = cent_mpc_fn(Nx, Nu, A, B, d, Q, S, tFIR, tSim, x0);

%% Calculate costs + plot 
obj    = get_cost_fn(Q, S, tSim, x, u);
objVal = get_cost_fn(Q, S, tSim, xVal, uVal);

% Output costs
fprintf('Distributed cost: %f\n', obj);
fprintf('Centralized cost: %f\n', objVal);

figure(1)
plot(1:tSim+1,xVal(1,:),'b',1:tSim+1,x(1,:),'*b',1:tSim+1,xVal(3,:),'g',1:tSim+1,x(3,:),'*g')
xlabel('$$Time$$','interpreter','latex','Fontsize', 10)
ylabel('$$\theta_{1},\ \theta_{2}$$','Interpreter','Latex','Fontsize', 10)
leg1 = legend('$$\theta_{1}\ Centralized\ MPC$$', '$$\theta_{1}\ Localized\ MPC\ using\ ADMM$$','$$\theta_{2}\ Centralized\ MPC$$', '$$\theta_{2}\ Localized\ MPC\ using\ ADMM$$');
set(leg1,'Interpreter','latex'); set(leg1, 'Fontsize', 8)
title('Subsystems 1 and 2')