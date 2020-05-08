clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User-defined
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model parameters
q = 0.8;
k = [2 4 0.1];
a = 1;

% initial conditions
x0 = [0.6; 0.6; 0.1];

% simulation length
tHorizon = 40;

% SLS length
tFIR = 20;

% disturbance (step)
val1 = 0;  % value of first disturbance
dur1 = 5;  % how long first disturbance lasts
val2 = -3; % value of second disturbance (lasts forever)
ws = [val1*ones(1, dur1) val2*ones(1, tHorizon-dur1)];

% sampling time for discretization
Ts = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup system
Nx = 3;
Nu = 3;

sys = LTISystem();
sys.Nx = Nx; sys.Nu = Nu;

xs      = zeros(Nx, tHorizon);
us      = zeros(Nu, tHorizon);
xs(:,1) = x0;

u_lb  = zeros(Nu, 1);
u_ub  = ones(Nu, 1);
x_lb  = [1e-3; 0.5; 1e-3];

for t=1:tHorizon-1
    fprintf('Time: %d\n', t);
    
    x_ = xs(:,t);
    u_ = us(:,t);
    w_ = ws(t);

    [Ac, Bc]        = linearize_glyco(x_, u_, q, k, a);    
    [sys.A, sys.B2] = discretize(Ac, Bc, Ts);
        
    % coordinate shift bounds
    u_tilde_lb  = u_lb - u_;
    u_tilde_ub  = u_ub - u_;
    x_tilde_lb  = x_lb - x_;
    
    y_shift = pinv(sys.A - eye(Nx)) * f_glyco(x_, u_, w_, q, k, a)*Ts;
    yt      = y_shift;
    y_lb    = x_tilde_lb + y_shift;
    
    relax = false;
    [y, u_tilde, status] = mpc_glyco(sys, tFIR, y_lb, u_tilde_lb, u_tilde_ub, yt, relax);

    % Due to linearization errors, sometimes a previous step
    % will yield a solution that does not satisfy the state bounds
    % in this case, try relaxing the constraint R{1}*xt >= x_lb 
    if strcmp(status, 'Infeasible')
        relax = true;
        [y, u_tilde, status] = mpc_glyco(sys, tFIR, y_lb, u_tilde_lb, u_tilde_ub, yt, relax);
    end
    
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

subplot(2,1,1);
hold on;
for i=1:2
    plot(1:tHorizon, xs(i,:));
end
ylabel('x1, x2');
legend('x1 (intermeds)', 'x2 (atp)');

subplot(2,1,2);
plot(1:tHorizon, xs(3,:));
ylabel('x3 (growth)');
