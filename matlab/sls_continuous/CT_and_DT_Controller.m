
% Use MOSEK

% Setup sls_master

%%
clc,clear,close all

setenv('PATH', [getenv('PATH') ';D:\Mosek\11.0\tools\platform\win64x86\bin']);
cvx_setup
cvx_solver mosek
% generate sys.A, sys.B2
% Nx = 3; alpha = 0.2; rho = 1; actDens = 1;
% sys = generate_dbl_stoch_chain(Nx, rho, actDens, alpha);


load('grid_model.mat');
[n,m] = size(sysB);
Q = eye(n);
% Q = diag([2,1,1]);
R = eye(m);
% R = diag([1,6,5]);


% 0.5*number of poles
hn = 3;

%% continuous-time model
Ac = full(sysA);
Bc = full(sysB);
Cc = eye(n);

ct_ol_sys = ss(Ac,Bc,Cc,0);

%% discrete-time model
dt_horizon = 2*hn+1;
Ts = 0.2;
dt_ol_sys = c2d(ct_ol_sys,0.1);
Ad = dt_ol_sys.A;
Bd = dt_ol_sys.B;
Cd = dt_ol_sys.C;

dtsys = LTISystem;
dtsys.A = Ad;
dtsys.B1 = eye(n);
dtsys.B2 = Bd;
dtsys.C1 = [Q^0.5;zeros(m,n)];
dtsys.D11 = zeros(m+n,n);
dtsys.D12 = [zeros(n,m);R^0.5];
dtsys.Nx = n;
dtsys.Nu = m;
dtsys.Nw = n;
dtsys.Nz = m+n;

dt_slsParams    = SLSParams();
dt_slsParams.T_ = dt_horizon;
dt_slsParams.add_objective(SLSObjective.H2, 1); % H2 objective, coefficient = 1

dt_clMaps  = state_fdbk_sls(dtsys, dt_slsParams);

dt_Phix = dt_clMaps.R_;
dt_Phiu = dt_clMaps.M_;

z = tf('z');
dt_Phix_tf = zeros(size(dt_Phix{1}));
dt_Phiu_tf = zeros(size(dt_Phiu{1}));

for ind = 1:dt_horizon
    dt_Phix_tf = dt_Phix_tf+dt_Phix{ind}/z^ind;
    dt_Phiu_tf = dt_Phiu_tf+dt_Phiu{ind}/z^ind;
end

dt_Phix1 = full(dt_Phix{1});
dt_Phix2 = full(dt_Phix{2});
dt_Phix3 = full(dt_Phix{3});
dt_Phix4 = full(dt_Phix{4});
dt_Phix5 = full(dt_Phix{5});
dt_Phix6 = full(dt_Phix{6});
dt_Phix7 = full(dt_Phix{7});

dt_Phiu1 = full(dt_Phiu{1});
dt_Phiu2 = full(dt_Phiu{2});
dt_Phiu3 = full(dt_Phiu{3});
dt_Phiu4 = full(dt_Phiu{4});
dt_Phiu5 = full(dt_Phiu{5});
dt_Phiu6 = full(dt_Phiu{6});
dt_Phiu7 = full(dt_Phiu{7});

err = norm((z*eye(n)-Ad)*dt_Phix_tf-Bd*dt_Phiu_tf-eye(n));
%% DT simulation
dt_on_ct_H2norm_sim = 0;
dt_H2norm_sim = 0;

Tsim = 1000;
dT = 0.002;
x0 = zeros(n,1);

for num = 1:n
x0 = zeros(n,1);
x0(num) = 1;
dt_sim_output = sim('dt_toy_model_simulation.slx');

dt_x = dt_sim_output.dt_x_sim.Data;
dt_u = squeeze(dt_sim_output.dt_u.Data)';

dt_on_ct_x = dt_sim_output.dt_on_ct_x_sim.Data;
dt_on_ct_u = squeeze(dt_sim_output.dt_on_ct_u.Data)';
dt_time = dt_sim_output.tout;

numberOfSteps = size(dt_x,1);
for ind = 1:numberOfSteps
    % Compute step length
    % if ind == 1
    %     steplength = ct_time(1);
    % else
    %     steplength = ct_time(ind)-ct_time(ind-1);
    % end
    steplength = dT;

    % add norm of the step
    dNorm = steplength*(dt_x(ind,:)*Q*dt_x(ind,:)'+dt_u(ind,:)*R*dt_u(ind,:)');
    dt_H2norm_sim = dt_H2norm_sim+dNorm;

    dNorm_on_ct = steplength*(dt_on_ct_x(ind,:)*Q*dt_on_ct_x(ind,:)'...
    +dt_on_ct_u(ind,:)*R*dt_on_ct_u(ind,:)');
    dt_on_ct_H2norm_sim = dt_on_ct_H2norm_sim+dNorm_on_ct;
end

end
dt_H2norm_sim = dt_H2norm_sim.^0.5
dt_on_ct_H2norm_sim = dt_on_ct_H2norm_sim.^0.5

%% Part II: LQR controller design. This is a "reference", our goal is to try and approach this value using SLS
[K_lqr,S_lqr,LQRPoles] = lqr(Ac,Bc,Q,R);
s = tf("s");
K_lqr = -K_lqr;
Phix_lqr = (s*eye(n)-Ac-Bc*K_lqr)^-1;
Phiu_lqr = K_lqr*Phix_lqr;
Psi_LQR = [Q^0.5 zeros(n,m); zeros(m,n), R^0.5] * [Phix_lqr; Phiu_lqr];

C_new = [Q^(1/2);R^(1/2)*K_lqr];
CLsys_lqr = ss(Ac+Bc*K_lqr,eye(n),C_new,0);
H2norm_LQR = norm(CLsys_lqr,2);
% LQRCost_2 = norm(Psi_LQR,2).^2;
disp("Based on LQR, the optimal H2 norm is:");
disp(H2norm_LQR);
disp('----------------------');
% disp('-----------');


%% Continuous-time SLS solution

hn = 3;

% This is for 4 poles:
ct_poles = getPoles(hn);

[ct_H2norm,ct_Psi,ct_Phix,ct_Phiu] = run_H2_cvx(Ac,Bc,Q,R,ct_poles,0,0);

ct_H2norm

% slsp1 = ct_poles(1);
% slsp2 = ct_poles(2);
% slsp3 = ct_poles(3);
% slsp4 = ct_poles(4);
% 
% Phix1 = ct_Phix{1};
% Phix2 = ct_Phix{2};
% Phix3 = ct_Phix{3};
% Phix4 = ct_Phix{4};
% 
% Phiu1 = ct_Phiu{1};
% Phiu2 = ct_Phiu{2};
% Phiu3 = ct_Phiu{3};
% Phiu4 = ct_Phiu{4};

%% Save data

% We want to save continuous-time model's Ac,Bc the weighting parameter Q,R,
% Ts (which determines what is Ad and Bd)

save('ct_and_dt_controller_parameters_1', 'Ac','Bc','Ad','Bd','Ts','hn', ...
    'dt_horizon','Q','R', 'H2norm_LQR',...
    'ct_Phix','ct_Phiu','dt_Phix','dt_Phiu','dt_on_ct_H2norm_sim','ct_H2norm');



%% Continuous-time H2 norm
% ct_H2norm_sim = 0;
% 
% Tsim = 60;
% dT = 0.002;
% 
% for num = 1:n
% x0 = zeros(n,1);
% x0(num) = 1;
% ct_sim_output = sim('ct_toy_model_simulation.slx');
% 
% ct_x = ct_sim_output.x_sim.Data;
% ct_u = squeeze(ct_sim_output.u.Data)';
% ct_time = ct_sim_output.tout;
% 
% numberOfSteps = size(ct_x,1);
% for ind = 1:numberOfSteps
%     % Compute step length
%     % if ind == 1
%     %     steplength = ct_time(1);
%     % else
%     %     steplength = ct_time(ind)-ct_time(ind-1);
%     % end
%     steplength = dT;
% 
%     % add norm of the step
%     dNorm = steplength*(ct_x(ind,:)*Q*ct_x(ind,:)'+ct_u(ind,:)*R*ct_u(ind,:)');
%     ct_H2norm_sim = ct_H2norm_sim+dNorm;
% end
% 
% end

%% Discrete-time solution

