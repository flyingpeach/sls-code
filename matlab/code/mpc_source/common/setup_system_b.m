% do not clear variables; takes numPendula, locality as input

%% general setup
rng(2020)
d = locality;
n = numPendula;

Nx = 2*n; 
Nu = n;

x0   = rand(Nx,1);
tSim = 10; % simulation time
tFIR = 5; % FIR time horizon

%% plant setup
Ac = zeros(Nx,Nx); j = 0; % A matrix (continuous)
for i = 1:2:Nx
    j = j+1;
    if j == 1
        Ac (i:i+1,i+2:i+3) = [0    0; 1  1];
        Ac (i:i+1,i:i+1) = [0 1; -3 -3];
    elseif j == Nx/2        
        Ac (i:i+1,i:i+1) = [0 1; -3 -3];
        Ac (i:i+1,i-2:i-1) = [0    0; 1 1];
    else
        Ac (i:i+1,i+2:i+3) = [0    0; 1 1];
        Ac (i:i+1,i:i+1) = [0 1; -3 -3];
        Ac (i:i+1,i-2:i-1) = [0    0; 1 1];
    end
end

% B matrix
Bc = zeros(Nx,Nu); j = 0;
for i = 1:2:Nx
    j = j+1;
    Bc (i:i+1,j) = [0; 1];
end

% Discretize 
Ts = .1;
A  = (eye(Nx)+Ac*Ts);
B  = Ts*Bc;