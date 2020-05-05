%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User-defined
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model parameter
q = 0.8;

% initial conditions
x0 = 0.1 * ones(Nx, 1);
u0 = zeros(Nu, 1);

% simulation length
tHorizon = 20;

% disturbance
ws = [0.1 * ones(1, tHorizon/2), ones(1, tHorizon/2)];

% sampling time for discretization
Ts = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup system
u_lb  = [-1; -1; -1];
u_ub  = [0; 1; 0];
x4_lb = 1e-3;
    
Nx       = 5;
Nu       = 3;

xs         = zeros(Nx, tHorizon);
us         = zeros(Nu, tHorizon);
xs(:,1) = x0; us(:,1) = u0;

for t=1:tHorizon-1
    x_ = xs(:,t);
    u_ = us(:,t); % always zero
    w_ = ws(:,t);

    [Ac, Bc] = linearize_bilo(x_, u_, q);
    [A, B]   = discretize(Ac, Bc, Ts);
    
    % coordinate shift bounds and pass in to MPC
    
    % do MPC with (A,B), x0=0 (since we're in coordinate-shifted regime)
    % remember to include x_(5) in objective    
    % MPC returns y (coordinate-shifted x), u_tilde 
    
    % calculate x(t+1) from dynamics and u_tilde directly
    x_tilde_nxt = f_bilo(x_, u_, w_, q)*Ts + B*u_tilde;
    xs(:,t+1) = x_ + x_tilde_next;
    
    % update actuation
    us(:,t) = u_ + u_tilde;    
end


