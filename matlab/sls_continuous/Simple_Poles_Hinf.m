% Use MOSEK

%% Part1: Generate Plant (This is not included in the first stage)
clc,clear,close all
% % Generate a chain system
% Nx = 4; 
% rho = 1;
% density = 1;
% alpha = 0.4;
% 
% % Two types of random chains
% % sys = generate_dbl_stoch_chain(Nx,rho,density,alpha);
% sys = generate_rand_chain(Nx,rho,density);
% 
% sysA = sys.A
% sysB = sys.B2


% Generate power grid system
% Properties of the Grid
sizeOfGrid = 3;
threshold = 0.8;
seed = 6;

% Generate matrices that shows properties of the system.
[AdjMtx, NodeCoords, SusceptMtx, InertiasInv, Dampings] = generate_grid_topology(sizeOfGrid,threshold,seed);
actuatedNodes = 1:sizeOfGrid^2;
powersys = generate_continuous_time_grid_plant(actuatedNodes,AdjMtx,SusceptMtx,InertiasInv,Dampings);
sysA = powersys.A;
sysB = powersys.B2;
plot_graph(AdjMtx,NodeCoords,'red');

save('small_grid.mat','sysA');
save('small_grid.mat','sysB','-append');

% % Load system
load('small_grid.mat');

% sysA = [1,0;0,1];
% sysB = [1,0;0,1];

% Set Q and R
n = size(sysB,1);
m = size(sysB,2);

Q = eye(n);
R = eye(m);

hnmin = 1;
hnmax = 3;

disp('A = ');
disp(sysA);
disp('B2 = ');
disp(sysB);
disp('----------------------');

%% Part II: hinfsyn controller design. This is a "reference", our goal is to try and approach this value using SLS

pA = full(sysA);
pB = [eye(n),full(sysB)];
pC = [Q^0.5;zeros(m,n);eye(n)];
pD = [zeros(n+m,n),[zeros(n,m);R^0.5];zeros(n),zeros(n,m)];

Plant = ss(pA,pB,pC,pD);

[Kctrl,CL,Hinf_ref] = hinfsyn(Plant,n,m);

Poles = eig(CL.A);
Hinf_ref;

disp("Based on Hinfsyn, the optimal Hinf norm is:");
disp(Hinf_ref);
disp('----------------------');

%% Part III: Verifying SPA (which is based on Archimedes pole approximation method)

SPA_SLSCost = [];
SPA_runtime = [];
SPA_statuses = [];
flag = 0; %Dense Controller

tic;
for hn = hnmin:hnmax
SLSPoles = getPoles(hn);
[Hinfnorm,Psi,Phix,Phiu,phix_supp,phiu_supp,status] = run_Hinf_cvx(sysA,sysB,Q,R,SLSPoles,flag,0);
newruntime=toc;
if(hn==hnmin)
    SPA_runtime = [SPA_runtime;newruntime];
else
SPA_runtime = [SPA_runtime;newruntime-SPA_runtime(size(SPA_runtime,1))];
end
SPA_SLSCost = [SPA_SLSCost;Hinfnorm];
SPA_statuses = [SPA_statuses,' '+status];

disp('Optimal solution is Obatined for hn being:');
disp(hn);
disp('Runtime for this iteration is')
disp(newruntime);
disp('- - - -');
end

disp('SPA Cost is:');
disp(SPA_SLSCost);
disp('----------------------');

% [flag,inds] = Hinf_sparsity_sanity_check(Phix,Phiu,phix_supp,phiu_supp)


%% Part IV: Effect of Sparsity

Sparse_SLSCost = [];
Sparse_runtime = [];
Sparse_statuses = [];

flag = 1; %grid model

hn = 2; %Fix set of poles
SLSPoles = getPoles(hn);

tic;
for d = 1:n
[Sparse_Hinf_norm,Psi,Phix,Phiu,phix_supp,phiu_supp,status] = run_Hinf_cvx(sysA,sysB,Q,R,SLSPoles,flag,d);
newruntime=toc;
Sparse_SLSCost = [Sparse_SLSCost;Sparse_Hinf_norm];
if(d==1)
    Sparse_runtime = [Sparse_runtime;newruntime];
else
Sparse_runtime = [Sparse_runtime;newruntime-Sparse_runtime(size(Sparse_runtime,1))];
end
Sparse_statuses = [Sparse_statuses,' '+status];
disp('Optimal solution is Obatined for d being:');
disp(d);
disp('Runtime for this iteration is')
disp(newruntime);
disp('- - - -');
end

disp('Cost Under different d is:');
disp(Sparse_SLSCost);
disp('----------------------');


%% Part V: Sparse SLS + SPA 
Sparse_SPA_SLSCost = [];
Sparse_SPA_runtime = [];
Sparse_SPA_statuses = [];

flag = 1;
d = 2; %Fix d = 1;

tic;
for hn = hnmin:hnmax
    SLSPoles = getPoles(hn);
    [Sparse_Hinf_SPA_norm,Psi,Phix,Phiu,phix_supp,phiu_supp,status] = run_Hinf_cvx(sysA,sysB,Q,R,SLSPoles,flag,d);
    newruntime=toc;
    Sparse_SPA_SLSCost = [Sparse_SPA_SLSCost;Sparse_Hinf_SPA_norm];
    if(hn==hnmin)
        Sparse_SPA_runtime = [Sparse_SPA_runtime;newruntime];
    else
        Sparse_SPA_runtime = [Sparse_SPA_runtime;newruntime-Sparse_SPA_runtime(size(Sparse_SPA_runtime,1))];
    end
    Sparse_SPA_statuses = [Sparse_SPA_statuses," "+status];
    disp('Optimal solution is Obatined for hn being:');
    disp(hn);
    disp('Runtime for this iteration is')
    disp(newruntime);
    disp('- - - -');
end

disp('----------------------');

disp('Sparse SPA Cost is:');
disp(Sparse_SPA_SLSCost);

disp('Runtime for SPA is:');
disp(SPA_runtime);

disp('Runtime under different d is');
disp(Sparse_runtime);

disp('Runtime for Sparse SPA is');
disp(Sparse_SPA_runtime);

%% Load the data
SPA_SLS_Ratio = SPA_SLSCost/Hinfnorm;
Sparse_SLS_Ratio = Sparse_SLSCost/Hinfnorm;
Sparse_SPA_Ratio = Sparse_SPA_SLSCost/Hinfnorm;
Dist = (1:n)';
PoleNumber = (2*hnmin:2:2*hnmax)';
