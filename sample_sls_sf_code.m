function sample_sls_sf_code();

%number of states
Nx = 10;
%FIR horizon
T = 20;
Tmax = 25;
Tsim = T + Tmax;

% Generate state space parameters
[A,B] = generate_dbl_stoch_chain(Nx,.2,1,1);

%number actuators
Nu = size(B,2);

% Specify objective function parameters
C = [speye(Nx); sparse(Nu,Nx)];
D = [sparse(Nx,Nu); speye(Nu)];

%% Solve for optimal system response (R,M) with H2 (squared) objective
% Basic means no constraints on (R,M) so this gives a centralized
% controller
[R,M]  = sf_sls_basic(A,B,C,D,T,'H2');

% Make a heatmap
%loc is where the disturbane hits
loc = floor(Nx/2);

make_heat_map(A,B,T,Nx,Nu,R,M,loc,Tsim,'Centralized')
% 
%% Solve for optimal system response (R,M) with H2 (squared) objective
% d-localized means d-hop locality constraints on (R,M)
% comm-speed needs to be sufficiently large.  ta is actuation delay
d = 3;
comms = 2;
ta = 1;

[R1,M1]  = sf_sls_d_localized(A,B,C,D,T,d,comms,ta,'H2');

% Make a heatmap
%loc is where the disturbane hits
loc = floor(Nx/2);
% Tmax is how long we run things for

make_heat_map(A,B,T,Nx,Nu,R1,M1,loc,Tsim,'Localized')

%% Solve for optimal system response (R,M) with H2 (squared) objective
% and approx d-localized means d-hop locality constraints on (R,M)
% robust_stab < 1 => we can guarantee 
d = 3;
comms = 1;
ta = 1;
lambda = 10^3;

[R2,M2,obj,robust_stab]  = sf_sls_approx_d_localized(A,B,C,D,T,d,comms,ta,lambda,'H2');

% Make a heatmap
%loc is where the disturbane hits
loc = floor(Nx/2);
% Tmax is how long we run things for

make_heat_map(A,B,T,Nx,Nu,R2,M2,loc,Tsim,'Approximately Localized')


