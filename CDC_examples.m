clc; clear; close all;

Nx = 50; % states
loc = floor(Nx/2); % specify where disturbance hits

actuation = 1; % actuation density
alpha = 0.4;
rho = 1;

% Construct system matrices
[A,B] = generate_dbl_stoch_chain(Nx,alpha,rho,actuation);
[~,Nu] = size(B); % number of actuators

% Specify objective function parameters
C = [speye(Nx); sparse(Nu,Nx)];
D = [sparse(Nx,Nu); speye(Nu)];

% locality
d = 6; 
comms = 2;
ta = 1;

%% Open Loop and Centralized Solutions
make_heat_map(A,B,10,Nx,Nu,0,0,loc,90,'open loop',1)

% Centralized Solution
T_centralized = 10;
Tmax = 25;
Tsim = T_centralized + Tmax;
[R_central,M_central]  = sf_sls_basic(A,B,C,D,T_centralized,'H2');
make_heat_map(A,B,T_centralized,Nx,Nu,R_central,M_central,loc,Tsim,[int2str(Nx), ' Node Centralized Solution'])

%% Vary Horizon Length
ta = 1;
T = [ 3,4,5,6,7,8]; % FIR horizon length
for i=1:length(T)
   Tmax = 16;
   Tsim = T(i) + Tmax;

   [R_T{i},M_T{i},obj_T(i)]  = sf_sls_d_localized(A,B,C,D,T(i),d,comms,ta,'H2');

   % Make a heatmap
   %loc is where the disturbane hits
   make_heat_map(A,B,T(i),Nx,Nu,R_T{i},M_T{i},loc,Tsim,[' T = ' int2str(T(i))])
end
figure
p = plot(T,obj_T,'o-');
title([int2str(Nx), ' Node Chain' ])
set ( gca, 'xdir', 'reverse' )
xlabel('FIR Horizon T')
ylabel('Localized H_2-Norm Cost')
set(gca,'FontSize',16,'fontWeight','bold')
set(p,'Color','red')
set(p,'LineWidth', 2)


%% Vary Actuator Density

actuation = [1 0.9 0.8 0.7 0.6 0.5 0.4 0.3];
ta = 1;
for i = 1:length(actuation)
   [A,B] = generate_dbl_stoch_chain(Nx,alpha,rho,actuation(i));
   [~,Nu] = size(B); % number of actuators 
   C = [speye(Nx); sparse(Nu,Nx)];
   D = [sparse(Nx,Nu); speye(Nu)];
   
   d = 6; 
   comms = 2;
   
   T = 15;
   Tmax = 45;
   Tsim = T + Tmax;
   [R_act{i},M_act{i},obj_Act(i)]  = sf_sls_d_localized(A,B,C,D,T,d,comms,ta,'H2');
   make_heat_map(A,B,T,Nx,Nu,R_act{i},M_act{i},loc,Tsim,['Actuation = ',num2str(actuation(i))])
end
figure
p1 = plot(flip(actuation),flip(obj_Act),'o-');
set ( gca, 'xdir', 'reverse' )
title([int2str(Nx), ' Node Chain' ])
xlabel('Actuation Density')
ylabel('Localized H_2-Norm Cost')
set(gca,'FontSize',16,'fontWeight','bold')
set(p1,'Color','red')
set(p1,'LineWidth', 2)

actuation(i+1) = 0.2; 
[A,B] = generate_dbl_stoch_chain(Nx,alpha,rho,actuation(i+1));
[~,Nu] = size(B); % number of actuators 
C = [speye(Nx); sparse(Nu,Nx)];
D = [sparse(Nx,Nu); speye(Nu)];
   
d = 6; 
comms = 2;

lambda = 10^3;
[R_rob_act,M_rob_act,obj_Act(i+1),robust_stab_act]  = sf_sls_approx_d_localized(A,B,C,D,T,d,comms,ta,lambda,'H2');
make_heat_map(A,B,T,Nx,Nu,R_rob_act,M_rob_act,loc,Tsim,['Actuation = ',num2str(actuation(i+1))])
figure

p1 = plot(flip(actuation),flip(obj_Act),'o-');
set ( gca, 'xdir', 'reverse' )
title([int2str(Nx), ' Node Chain' ])
xlabel('Actuation Density')
ylabel('Localized H_2-Norm Cost')
set(gca,'FontSize',16,'fontWeight','bold')
set(p1,'Color','red')
set(p1,'LineWidth', 2)

%% Vary d

ta = 1;
T = 10;
d = [10 8 6 5 4 3];
Tmax = 16;
Tsim = T + Tmax;
actuation = 1;
[A,B] = generate_dbl_stoch_chain(Nx,alpha,rho,actuation);
[~,Nu] = size(B); % number of actuators
C = [speye(Nx); sparse(Nu,Nx)];
D = [sparse(Nx,Nu); speye(Nu)];

for i=1:length(d)
   [R_d{i},M_d{i},obj_d(i)]  = sf_sls_d_localized(A,B,C,D,T,d(i),comms,ta,'H2');

   % Make a heatmap
   %loc is where the disturbane hits
   make_heat_map(A,B,T,Nx,Nu,R_d{i},M_d{i},loc,Tsim,['d = ',int2str(d(i))])
end

d(i+1) = 2;
[R_rob_d,M_rob_d,obj_d(i+1),robust_stab_d]  = sf_sls_approx_d_localized(A,B,C,D,T,d(i+1),comms,ta,lambda,'H2');
make_heat_map(A,B,T,Nx,Nu,R_rob_d,M_rob_d,loc,Tsim,['d = ',int2str(d(i+1))])

figure
p2 = plot(d,obj_d,'o-');
set ( gca, 'xdir', 'reverse' )
title([int2str(Nx), ' Node Chain' ])
xlabel('d-hops')
ylabel('Localized H_2-Norm Cost')
set(gca,'FontSize',16,'fontWeight','bold')
set(p2,'Color','red')
set(p2,'LineWidth', 2)

%% Comm speed

 

T = 10;
d = 8;
comms = [4 3 2 1.75 1.5 1.25 1];
for i=1:length(comms)
   Tmax = 25;
   Tsim = T + Tmax;

   [R_comms{i},M_comms{i},obj_comms(i)]  = sf_sls_d_localized(A,B,C,D,T,d,comms(i),ta,'H2');

   % Make a heatmap
   %loc is where the disturbane hits
   make_heat_map(A,B,T,Nx,Nu,R_comms{i},M_comms{i},loc,Tsim,[' Comms = ' int2str(comms(i))])
end

figure
p4 = plot(comms,obj_comms,'o-');
set ( gca, 'xdir', 'reverse' )
title([int2str(Nx), ' Node Chain' ])
xlabel('\alpha')
ylabel('Localized H_2-Norm Cost')
set(gca,'FontSize',16,'fontWeight','bold')
set(p4,'Color','red')
set(p4,'LineWidth', 2)




%% Spread of state
d = 6;
comms = 2;
alpha = linspace(0,0.8,10);
actuation = 1; % actuation density
rho = 1;
T = 10;

for i=1:length(alpha)
   % Construct system matrices
   [A,B] = generate_dbl_stoch_chain(Nx,alpha(i),rho,actuation);
   disp(['Alpha = ', num2str(alpha(i))])
   disp((['Spectral Radius:']))
   max(abs(eig(full(A))))
   [~,Nu] = size(B); % number of actuators

   % Specify objective function parameters
   C = [speye(Nx); sparse(Nu,Nx)];
   D = [sparse(Nx,Nu); speye(Nu)];
   
    [R_alpha{i},M_alpha{i},obj_alpha(i)]  = sf_sls_d_localized(A,B,C,D,T,d,comms,ta,'H2');
    make_heat_map(A,B,T,Nx,Nu,R_alpha{i},M_alpha{i},loc,Tsim,['\alpha = ' num2str(alpha(i))])
end

figure
p3 = plot(alpha,obj_alpha,'o-');
title([int2str(Nx), ' Node Chain' ])
xlabel('\alpha')
ylabel('Localized H_2-Norm Cost')
set(gca,'FontSize',16,'fontWeight','bold')
set(p3,'Color','red')
set(p3,'LineWidth', 2)