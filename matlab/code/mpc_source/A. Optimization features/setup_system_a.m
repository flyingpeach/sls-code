clear all; close all; clc;

% general setup
rng(2020)
d = 3; % locality (this is between states, not subsystems)
n = 4; % number of pendulums

Nx = 2*n; 
Nu = n;

x0   = rand(Nx,1);
tSim = 30; % simulation time
tFIR = 20; % FIR time horizon

% plant setup
m = 1; k = 1; g = 10; l = 1;
block_off_diag  = [0    0; k*l/m  d/(m*l)];
block_diag_extr = [0 1; -g-k*l/m -d/(m*l)];
block_diag      = [0 1; -g-2*k*l/m -2*d/(m*l)];

Ac = zeros(Nx,Nx); j = 0; % A matrix (continuous time)
for i = 1:2:Nx
    j = j+1;
    if j == 1 % first node
        Ac (i:i+1,i+2:i+3) = block_off_diag;
        Ac (i:i+1,i:i+1) = block_diag;
    elseif j == Nx/2  % last node      
        Ac (i:i+1,i:i+1) = block_diag;
        Ac (i:i+1,i-2:i-1) = block_off_diag;
    else
        Ac (i:i+1,i+2:i+3) = block_off_diag;
        Ac (i:i+1,i:i+1) = block_diag;
        Ac (i:i+1,i-2:i-1) = block_off_diag;
    end
end

% B matrix (continous time)
Bc = zeros(Nx, Nu); j = 0;
for i = 1:2:Nx
    j = j+1;
    Bc (i:i+1,j) = [0; 1];
end

% Discretize 
Ts = .1;
A  = (eye(Nx)+Ac*Ts);
B  = Ts*Bc;