clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User-defined
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model parameters
q = 0.8;
k = 1;
a = 1;
h = 3;
g = 0.1;

% initial conditions
x0 = [1; 1];

% simulation length
tHorizon = 400;

% SLS length
tFIR = 4;

% disturbance (step)
ws = -1.1 * ones(1, tHorizon);

% sampling time for discretization
Ts = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup
Nx = 2;
Nu = 2;

sys = LTISystem();
sys.Nx = Nx; sys.Nu = Nu;

xs      = zeros(Nx, tHorizon);
us      = zeros(Nu, tHorizon);
xs(:,1) = x0;

u_ub = zeros(Nu, 1);
x_lb = 1e-3 * ones(Nx, 1);

for t=1:tHorizon-1
    fprintf('Time: %d\n', t);

    x=xs(:,t); x1=x(1); x2=x(2);
    w=ws(t);
        
    u1 = -log(1 + x2.^(2*h));
    u2 = -log(1 + x2.^(2*g));
    u  = [u1; u2];
    
    xs(:, t+1) = x + f_glyco(x, u, w, q, k, a) * Ts ;    
end

%% Plot
figure(); 
hold on;

for i=1:2
    plot(1:tHorizon, xs(i,:));
end

legend('x1 (intermeds)', 'x2 (atp)');