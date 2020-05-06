clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User-defined
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model parameters
q = 0.8;
k = 1;
a = 1;

% initial conditions
x0 = [1.1; 1.1];

% simulation length
tHorizon = 30;

% SLS length
tFIR = 8;

% disturbance (step)
ws = -1.1 * ones(1, tHorizon);

% sampling time for discretization
Ts = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup system
Nx = 2;
Nu = 2;

sys = LTISystem();
sys.Nx = Nx; sys.Nu = Nu;

xs      = zeros(Nx, tHorizon);
us      = zeros(Nu, tHorizon);
xs(:,1) = x0;

u_ub  = zeros(Nu, 1);
x_lb  = 1e-3 * ones(Nx, 1);
x_ref = ones(Nx,1);

for t=1:tHorizon-1
    fprintf('Time: %d\n', t);
    
    x_ = xs(:,t);
    u_ = us(:,t);
    w_ = ws(t);

    [Ac, Bc]        = linearize_glyco(x_, u_, q, k, a);    
    [sys.A, sys.B2] = discretize(Ac, Bc, Ts);
    
    % coordinate shift bounds
    u_tilde_ub  = u_ub - u_;
    x_tilde_lb  = x_lb - x_;
    x_tilde_ref = x_ref - x_;
    
    y_shift = pinv(sys.A - eye(Nx)) * f_glyco(x_, u_, w_, q, k, a)*Ts;
    yt      = y_shift;
    y_lb    = x_tilde_lb + y_shift;
    y_ref   = x_tilde_ref + y_shift;
    
    [y, u_tilde] = mpc_glyco(sys, tFIR, y_lb, u_tilde_ub, yt, y_ref);
    
    % calculate x(t+1) from dynamics and u_tilde directly
    x_tilde_nxt = f_glyco(x_, u_, w_, q, k, a)*Ts + sys.B2*u_tilde;
    xs(:,t+1)   = x_ + x_tilde_nxt;
    
    % update actuation
    us(:,t)   = u_ + u_tilde;
    % next time, linearize about current u (more accurate than zeros)
    us(:,t+1) = us(:,t);

end

%% Plot
figure(); 
hold on;

for i=1:2
    plot(1:tHorizon, xs(i,:));
end

legend('x1 (intermeds)', 'x2 (atp)');
