% Algorithm II, System A

%% Setup
setup_system_a;

% Weights
Q = zeros(Nx,Nx);
for i = 1:2:Nx
    Q(i,i) = 1;
    if i > 1
        Q(i,i-2) = -1/2;
    end
    if i < Nx-1
        Q(i,i+2) = -1/2;
    end
    if i < Nx
        Q(i+1,i+1) = .01;
    end
end
S = eye(Nu);

% Build cost matrix
C = [];
for t = 0:tFIR-1
    C = blkdiag(C, Q);
end
for t = 0:tFIR-2
    C = blkdiag(C, S);
end    

% Coupling
indices = cell(1, length(C));
for i = 1:length(C)
    for j = 1:length(C)
        if C(i,j) ~= 0
            indices{i} = [indices{i} j];
        end
    end
end

%% Distributed MPC
rho   = 1; 
eps_d = 1e-3; % convergence criterion for ||Psi(k+1) - Psi(k)||
eps_p = 1e-4; % convergence criterion for ||Phi(k+1) - Psi(k+1)||
maxIters = 5000;

mu       = 1;
eps_pres = 1e-3;
eps_dres = 1e-4;
maxConsensusIters = 200;

[x, u, ~] = mpc_algorithm2(Nx, Nu, A, B, d, tFIR, tSim, x0, ... % system
                              eps_d, eps_p, rho, maxIters, ...
                              indices, eps_pres, eps_dres, mu, ...
                              maxConsensusIters, C);

 %% Centralized MPC (for validation + comparison)
[xVal, uVal, ~] = mpc_centralized(Nx, Nu, A, B, d, Q, S, tFIR, tSim, x0);

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