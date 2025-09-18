% This file includes all the process that are of interest

% MOSEK or Gurobi, Gurobi is probably faster


%% Part1: Generate Plant (This is not included in the first stage)  
clc,clear,close all
% sysA = [1,0.1,0;0.2,1,0.2;0,-0.3,1]
% sysB = eye(3)
% % Generate a chain system
% Nx = 5; 
% rho = 1;
% density = 1;
% alpha = 0.4;
% 
% % Two types of random chains
% sys = generate_dbl_stoch_chain(Nx,rho,density,alpha);
% % sys = generate_rand_chain(Nx,rho,density);
% 
% sysA = sys.A;
% sysB = sys.B2;

% % Generate power grid system
% % Properties of the Grid
% sizeOfGrid = 3;
% threshold = 0.7;
% seed = 18;
% 
% % Generate matrices that shows properties of the system.
% [AdjMtx, NodeCoords, SusceptMtx, InertiasInv, Dampings] = generate_grid_topology(sizeOfGrid,threshold,seed);
% actuatedNodes = 1:sizeOfGrid^2;
% powersys = generate_continuous_time_grid_plant(actuatedNodes,AdjMtx,SusceptMtx,InertiasInv,Dampings);
% plot_graph(AdjMtx,NodeCoords,'red');
% sysA = powersys.A;
% sysB = powersys.B2;
% 
% save('grid_model.mat','sysA');
% save('grid_model.mat','sysB','-append');

% sysA = [1,2,0;3,-5,6;9,2,3];
% sysB = [1;0;0];
% sysA = [1,2;3,4];
% sysB = [2;-1];
% imagesc(sysB ~=0)

% sysA = [1,0;0,1];
% sysB = [1,0;0,1];


% Load system
load('grid_model.mat');

SparseMode = 1; %Grid system - 1, Chain system - 2, Dense - 0.

% Set Q and R
n = size(sysB,1);
m = size(sysB,2);

Q = eye(n);
R = eye(m);

hnmin = 1;
hnmax = 4;


disp('A = ');
disp(sysA);
disp('B2 = ');
disp(sysB);
disp('----------------------');

% [~,~,Phix,Phiu,~,~,~] = run_H2_cvx(sysA,sysB,Q,R,-2,SparseMode,1)

%% Part II: LQR controller design. This is a "reference", our goal is to try and approach this value using SLS
[K_lqr,S_lqr,LQRPoles] = lqr(full(sysA),full(sysB),Q,R);
s = tf("s");
K_lqr = -K_lqr;
Phix_lqr = (s*eye(n)-sysA-sysB*K_lqr)^-1;
Phiu_lqr = K_lqr*Phix_lqr;
Psi_LQR = [Q^0.5 zeros(n,m); zeros(m,n), R^0.5] * [Phix_lqr; Phiu_lqr];

C_new = [Q^(1/2);R^(1/2)*K_lqr];
CLsys_lqr = ss(sysA+sysB*K_lqr,eye(n),C_new,0);
H2norm_LQR = norm(CLsys_lqr,2);
% LQRCost_2 = norm(Psi_LQR,2).^2 
disp("Based on LQR, the optimal H2 norm is:");
disp(H2norm_LQR);
disp('----------------------');
% disp('-----------');

%% Part III: Verifying SPA (which is based on Archimedes pole approximation method)
SPA_SLSCost = [];
SPA_runtime = [];
SPA_statuses = [];
flag = 0; %Dense Controller

tic;
for hn = hnmin:hnmax
SLSPoles = getPoles(hn);
% SLSPoles = -1:-1:-hn;
[H2_norm,Psi,Phix,Phiu,phix_supp,phiu_supp,status] = run_H2_cvx(sysA,sysB,Q,R,SLSPoles,flag,0);
newruntime=toc;
if(hn==hnmin)
    SPA_runtime = [SPA_runtime;newruntime];
else
SPA_runtime = [SPA_runtime;newruntime-SPA_runtime(size(SPA_runtime,1))];
end
SPA_SLSCost = [SPA_SLSCost;H2_norm];
SPA_statuses = [SPA_statuses;status];
disp('Optimal solution is Obatined for hn being:');
disp(hn);
disp('Runtime for this iteration is')
disp(newruntime);
disp('- - - -');
end

disp('SPA Cost is:');
disp(SPA_SLSCost);
disp('----------------------');

%% Part IV: Effect of Sparsity

Sparse_SLSCost = [];
Sparse_runtime = [];
Sparse_statuses = [];

flag = SparseMode; %grid model

hn = 3; %Fix set of poles
SLSPoles = getPoles(hn);

tic;
Dist = 1:n/2;
for d = 1:n/2
[H2_norm,Psi,Phix,Phiu,phix_supp,phiu_supp,status] = run_H2_cvx(sysA,sysB,Q,R,SLSPoles,flag,d);
newruntime=toc;
Sparse_SLSCost = [Sparse_SLSCost;H2_norm];
if(d==1)
    Sparse_runtime = [Sparse_runtime;newruntime];
else
Sparse_runtime = [Sparse_runtime;newruntime-Sparse_runtime(size(Sparse_runtime,1))];
end
Sparse_statuses = [Sparse_statuses;status];

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

flag = SparseMode;
d = 2; %Fix d 
PoleNumber = [2*hnmin:2:2*hnmax];

tic;
for hn = hnmin:hnmax
    SLSPoles = getPoles(hn);
    [H2_norm,Psi,Phix,Phiu,phix_supp,phiu_supp,status] = run_H2_cvx(sysA,sysB,Q,R,SLSPoles,flag,d);
    newruntime=toc;
    Sparse_SPA_SLSCost = [Sparse_SPA_SLSCost;H2_norm];
    if(hn==hnmin)
        Sparse_SPA_runtime = [Sparse_SPA_runtime;newruntime];
    else
        Sparse_SPA_runtime = [Sparse_SPA_runtime;newruntime-Sparse_SPA_runtime(size(Sparse_SPA_runtime,1))];
    end
    Sparse_SPA_statuses = [Sparse_SPA_statuses;status];
    disp('Optimal solution is Obatined for hn being:');
    disp(hn);
    disp('Runtime for this iteration is')
    disp(newruntime);
    disp('- - - -');
end



%%
disp('----------------------');
disp('Sparse SPA Cost is:');
disp(Sparse_SPA_SLSCost);

disp('----------------------');
disp('Runtime for SPA is:');
disp(SPA_runtime);

disp('Runtime under different d is');
disp(Sparse_runtime);

disp('Runtime for Sparse SPA is');
disp(Sparse_SPA_runtime);

%% Load Data
SPA_SLS_Ratio = SPA_SLSCost/H2norm_LQR;
Sparse_SLS_Ratio = Sparse_SLSCost/H2norm_LQR;
Sparse_SPA_Ratio = Sparse_SPA_SLSCost/H2norm_LQR;

