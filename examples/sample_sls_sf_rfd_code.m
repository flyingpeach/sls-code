function sample_sls_sf_rfd_code();

%number of states
Nx = 10;
%FIR horizon
T = 15;
Tmax = 25;
Tsim = T + Tmax;

randn('seed',0);
% Generate state space parameters
[A,B] = generate_rand_chain(Nx,.8,1);

%number actuators
Nu = size(B,2);

% Specify objective function parameters
C = [speye(Nx); sparse(Nu,Nx)];
D = [sparse(Nx,Nu); speye(Nu)];

%% Centralized Actuator RFD with H2 performance
num_acts = [];
clnorms = [];
for power = -2:1:3
    lambda = 10^power;
    [R,M,acts]  = sf_sls_basic_rfd(A,B,C,D,T,lambda,'H2');    
    Bd = B(:,acts);
    Dd = D(:,acts);
    %polishing_step
    [R,M,clnorm] = sf_sls_basic(A,Bd,C,Dd,T,'H2');
    num_acts = [num_acts; length(acts)];
    clnorms = [clnorms; clnorm];
end
num_acts
clnorms
figure
plot(num_acts,clnorms,'*-');
xlabel('Number of actuators')
ylabel('Close loop norm')
title('Centralized RFD tradeoff curve')

%% d-localized Actuator RFD with H2 performance
d = 3;
comms = 2;
ta = 1;

num_acts = [];
clnorms = [];
for power = -2:1:3
    lambda = 10^power;
    [R,M,acts]  = sf_sls_d_localized_rfd(A,B,C,D,T,d,comms,ta,lambda,'H2');    
    Bd = B(:,acts);
    Dd = D(:,acts);
    %polishing_step
    [R,M,clnorm]  = sf_sls_d_localized(A,Bd,C,Dd,T,d,comms,ta,'H2');
    num_acts = [num_acts; length(acts)];
    clnorms = [clnorms; clnorm];
end
num_acts
clnorms
figure
plot(num_acts,clnorms);
plot(num_acts,clnorms,'*-');
xlabel('Number of actuators')
ylabel('Close loop norm')
title('d-localized RFD tradeoff curve')


%% approx d-localized Actuator RFD with H2 performance
d = 3;
comms = 2;
ta = 1;

num_acts = [];
clnorms = [];
lambda_delta = 10^4;
for power = -2:1:3
    lambda = 10^power;
    [R,M,acts]  = sf_sls_approx_d_localized_rfd(A,B,C,D,T,d,comms,ta,lambda,lambda_delta,'H2');    
    Bd = B(:,acts);
    Dd = D(:,acts);
    %polishing_step
    [R,M,clnorm]  = sf_sls_approx_d_localized(A,Bd,C,Dd,T,d,comms,ta,lambda_delta,'H2');
    num_acts = [num_acts; length(acts)];
    clnorms = [clnorms; clnorm];
end
num_acts
clnorms
figure
plot(num_acts,clnorms);
plot(num_acts,clnorms,'*-');
xlabel('Number of actuators')
ylabel('Close loop norm')
title('Approx d-localized RFD tradeoff curve')



