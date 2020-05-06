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
tHorizon = 200;

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

for t=1:tHorizon-1
    fprintf('Time: %d\n', t);
    
    x_ = xs(:,t);
    u_ = us(:,t);
    w_ = ws(t);

    [Ac, Bc]        = linearize_glyco(x_, u_, q, k, a);    
    [sys.A, sys.B2] = discretize(Ac, Bc, Ts);

    x1=x_(1); x2=x_(2);
    u1 = -log(1 + x2.^(2*h));
    u2 = -log(1 + x2.^(2*g));
    u  = [u1; u2];

    u_tilde = u - u_;
    
    % calculate x(t+1) from dynamics and u_tilde directly
    x_tilde_nxt = f_glyco(x_, u, w_, q, k, a)*Ts+ sys.B2*u_tilde;
    xs(:,t+1)   = x_ + x_tilde_nxt;
    
    % update actuation
    us(:,t) = u;    
    % next time, linearize about current u (more accurate than zeros)
    us(:,t+1) = u;
end

%% Plot
figure(); 
hold on;

for i=1:2
    plot(1:tHorizon, xs(i,:));
end

legend('x1 (intermeds)', 'x2 (atp)');