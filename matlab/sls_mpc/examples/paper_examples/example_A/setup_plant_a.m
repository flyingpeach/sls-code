function sys = setup_plant_a()

n  = 4; % number of pendulums
Nx = 2*n; 
Nu = n;

m = 1; k = 1; f = 3; g = 10; l = 1;
block_off_diag  = [0    0; k*l/m  f/(m*l)];
block_diag      = [0 1; -g-2*k*l/m -2*f/(m*l)];

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

% Discretize + set up system
Ts = .1;

sys     = LTISystem();
sys.Nx  = Nx;
sys.Nu  = Nu;
sys.A   = (eye(Nx)+Ac*Ts);
sys.B2  = Ts*Bc;

sys.sanity_check_mpc();
end