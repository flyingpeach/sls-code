function sys = setup_plant_b(numPendula)

n  = numPendula;
Nx = 2*n; 
Nu = n;

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

% B matrix (continuous time)
Bc = zeros(Nx,Nu); j = 0;
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