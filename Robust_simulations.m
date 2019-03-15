clc; clear; close all;

Nx = 50; % states
loc = floor(Nx/2); % specify where disturbance hits

actuation = 0.5; % actuation density
alpha = 0.4;
rho = 1;

% Construct system matrices
[A,B] = generate_dbl_stoch_chain(Nx,alpha,rho,actuation);
[~,Nu] = size(B); % number of actuators

% Specify objective function parameters
C = [speye(Nx); sparse(Nu,Nx)];
D = [sparse(Nx,Nu); speye(Nu)];

T = 10;

% locality
d = 6; 
ta = 1;

Tmax = 25;
Tsim = T + Tmax;

comms = [2 1.5 1.4 1.3 1.2 1.1 1 0.9 0.8 0.7 0.6 0.5 0.4];
lambda = 1000;
for i=1:length(comms)
   [R2,M2,obj(i),robust_stab(i)]  = sf_sls_approx_d_localized(A,B,C,D,T,d,comms(i),ta,lambda,'H2');
    make_heat_map(A,B,T,Nx,Nu,R2,M2,loc,Tsim,['Comms = ',num2str(comms(i))])
end

figure
p1 = plot(comms,obj,'o-');
set ( gca, 'xdir', 'reverse' )
title([int2str(Nx), ' Node Chain' ])
xlabel('Comms Speed')
ylabel('Localized H_2-Norm Cost')
set(gca,'FontSize',16,'fontWeight','bold')
set(p1,'Color','red')
set(p1,'LineWidth', 2)

figure
p2 = plot(comms,robust_stab,'o-');
set ( gca, 'xdir', 'reverse' )
title([int2str(Nx), ' Node Chain' ])
xlabel('Comms Speed')
ylabel('Stability Margin')
set(gca,'FontSize',16,'fontWeight','bold')
set(p2,'Color','red')
set(p2,'LineWidth', 2)

