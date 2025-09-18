% Use Mosek

%% Part1: Generate Plant (This is not included in the first stage)  
clc,clear,close all
% Generate a chain system
Nx = 11; 
rho = 0.6;
density = 1;
alpha = 2/3;

sys = generate_dbl_stoch_chain(Nx,rho,density,alpha);

sysA = sys.A;
sysB = sys.B2;

SparseMode = 2; %Chain system

% Set Q and R
n = size(sysB,1);
m = size(sysB,2);

Q = eye(n);
R = eye(m);

% Display A and B matrix for the system
disp('A = ');
disp(sysA);
disp('B2 = ');
disp(sysB);
disp('----------------------');

%% Part II: LQR controller design. This is a "reference", our goal is to try and approach this value using SLS
[K_lqr,S_lqr,LQRPoles] = lqr(full(sysA),full(sysB),Q,R);
s = tf("s");
K_lqr = -K_lqr;
Phix_lqr = (s*eye(n)-sysA-sysB*K_lqr)^-1;
Phiu_lqr = K_lqr*Phix_lqr;

C_new = [Q^(1/2);R^(1/2)*K_lqr];
CLsys_lqr = ss(sysA+sysB*K_lqr,eye(n),C_new,0);
H2norm_LQR = norm(CLsys_lqr,2);
disp("Based on Hinfsyn, the optimal Hinf norm is:");
disp(H2norm_LQR);
disp('----------------------');

%% SLS Controllers

d = 2; % Communication distance
hn = 2; % Pole number = hn*2
SLSPoles = getPoles(hn);
[H2_norm,Psi,Phix,Phiu,phix_supp,phiu_supp,status] = run_H2_cvx(sysA,sysB,Q,R,SLSPoles,SparseMode,d);
H2_norm = H2_norm.^0.5;

% These variables are for simulation
Phix1 = full(Phix{1});
Phix2 = full(Phix{2});
Phix3 = full(Phix{3});
Phix4 = full(Phix{4});

Phiu1 = full(Phiu{1});
Phiu2 = full(Phiu{2});
Phiu3 = full(Phiu{3});
Phiu4 = full(Phiu{4});

slsp1 = SLSPoles(1);
slsp2 = SLSPoles(2);
slsp3 = SLSPoles(3);
slsp4 = SLSPoles(4);

sim('SLS_ChainModel_Sim');

Xsim = ans.x_sim.Data;
BUsim = squeeze(ans.u.Data);

Tsim = ans.tout;

ratio = H2_norm/H2norm_LQR


%% Plotting the heat map

Phix_sp = ((Phix1~=0)|(Phix2~=0)|(Phix3~=0)|(Phix4~=0));
Phiu_sp = ((Phiu1~=0)|(Phiu2~=0)|(Phiu3~=0)|(Phiu4~=0));

f2 = figure(2);
lmin = 0;
lmax = 1;
subplot(1,2,1)
clim([lmin,lmax]);
imagesc(Phix_sp)
title('$\mathbf{\Phi}_x(s)$', 'Interpreter', 'latex', 'FontSize', 20)
set(gca, 'XTick', [], 'YTick', []);
pos = get(gca, 'Position');
set(gca, 'Position', [pos(1)-0.07 pos(2) pos(3)+0.06 pos(4)]);
subplot(1,2,2)
imagesc(Phiu_sp);

title('$\mathbf{\Phi}_u(s)$', 'Interpreter', 'latex', 'FontSize', 20)
pos = get(gca, 'Position');
set(gca, 'Position', [pos(1)-0.02 pos(2) pos(3)+0.06 pos(4)]);
set(gca, 'XTick', [], 'YTick', []);
set(f2,'position',[100,100,700,345])

plot_ct_heat_map(Xsim',BUsim,'');

function plot_ct_heat_map(x, Bu, myTitle)
% Plots log-based heat map for x, u
% Inputs
%    x, Bu    : state and actuation values at nodes
%    myTitle  : overall title of the heat maps

f1 = figure;
logmin = -5;
logmax = 0;
fs = 15;

n = 256;                        % number of colors
c1 = [1, 1, 1];                 % white (low end)
c2 = [0, .1, .5];                 % blue (high end)
cmap = [linspace(c1(1),c2(1),n)', linspace(c1(2),c2(2),n)', linspace(c1(3),c2(3),n)'];
colormap(cmap);

subplot(1,2,1)
imagesc(log10(abs(x))); colorbar; colormap(cmap); clim([logmin logmax]);
title([myTitle ' $\mathbf{log_{10}|x|}$'], 'Interpreter', 'latex', 'FontSize',20); 
xlabel('\textbf{Time (s)}', 'Interpreter', 'latex', 'FontSize',18); 
ylabel('\textbf{Space}', 'Interpreter', 'latex', 'FontSize',18);
%pos = get(gca, 'Position');
%set(gca, 'Position', [pos(1)-0.05 pos(2)+0.02 pos(3)-0.00 pos(4)-0.05]);
ax = gca;
xTicks = ax.XTick;
ax.XTickLabel = xTicks/1000;
ax.XAxis.FontSize = fs;
ax.YAxis.FontSize = fs;

subplot(1,2,2)
imagesc(log10(abs(Bu))); colorbar; colormap(cmap); clim([logmin logmax]);
c = colorbar; c.FontSize = fs;
title([myTitle ' $\mathbf{log_{10}|u|}$'], 'Interpreter', 'latex', 'FontSize',20); 
xlabel('\textbf{Time (s)}', 'Interpreter', 'latex', 'FontSize',18); 
%pos = get(gca, 'Position');
%set(gca, 'Position', [pos(1)-0.02 pos(2)+0.02 pos(3)+0.00 pos(4)-0.05]);
ax = gca;
xTicks = ax.XTick;
ax.XTickLabel = xTicks/1000;
ax.XAxis.FontSize = fs;
ax.YAxis.FontSize = fs;

% These are the font settings from previous papers; uncomment as wanted
% set(gca,'FontSize',16,'fontWeight','bold')
% set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','bold')
set(f1,'position',[100,100,800,350])
end

