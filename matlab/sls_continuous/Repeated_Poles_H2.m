clc,clear,close all

% Use MOSEK

% Before running this code, update the saved path at the end.

%% SLS related parameters: pole location and maximum multiplicity
p = -0.3;
mult = 3;
% Q and R can be selected once A and B are determined

%% Define system, control objective (Q and R) and SLS parameters

% % For simple model, refer to these codes.
% A = 1;
% B = 2;
% 
% A = [1,3;-4,5];
% B = [1,0;0,1];
% A = [1,0.1,0;0.2,1,0.2;0,-0.3,1];
% B = eye(3);

% % For chain model, refer to these codes.
% Nx = 5; 
% rho = 1;
% density = 1;
% alpha = 0.4;
% 
% % Two types of random chains
% sys_chain = generate_dbl_stoch_chain(Nx,rho,density,alpha);
% % sys = generate_rand_chain(Nx,rho,density);
% 
% A = sys_chain.A;
% B = sys_chain.B2;

% For grid model, refer to these codes
load('grid_model.mat')

A = full(sysA);
B = full(sysB);

% System parameters
[n,m] = size(B);

Q = eye(n);
R = eye(m);



%% Part II: LQR controller design. This is a "reference", our goal is to try and approach this value using SLS
[K_lqr,S_lqr,LQRPoles] = lqr(A,B,Q,R);
s = tf("s");
K_lqr = -K_lqr;
Phix_lqr = (s*eye(n)-A-B*K_lqr)^-1;
Phiu_lqr = K_lqr*Phix_lqr;
Psi_LQR = [Q^0.5 zeros(n,m); zeros(m,n), R^0.5] * [Phix_lqr; Phiu_lqr];

C_new = [Q^(1/2);R^(1/2)*K_lqr];
CLsys_lqr = ss(A+B*K_lqr,eye(n),C_new,0);
H2norm_LQR = norm(CLsys_lqr,2);
% LQRCost_2 = norm(Psi_LQR,2).^2;
disp("Based on LQR, the optimal H2 norm is:");
disp(H2norm_LQR);
disp('----------------------');
% disp('-----------');

%% Solving for optimal controller using cvx
% Forming Ar and Br based on A,B and maximum multiplicity
Ar = zeros(mult);
Br = ones(mult,1);

for ind = 1:mult
    if ind<mult
        Ar(1+(ind-1)*1:ind*1,1+(ind-1)*1:(1+ind)*1) = [p,1];
    else
        Ar(1+(ind-1)*1:ind*1,1+(ind-1)*1:(ind)*1) = p*eye(1);
    end
end

Phix = cell(mult,1);
Phiu = cell(mult,1);
Psi = cell(mult,1);
Cr = cell(n+m,n);
Pr = cell(n+m,n);

cvx_begin sdp

variable phix(n,n,mult)
variable phiu(m,n,mult)
expression C_stack(m+n,n,mult)

variable P(m+n,n,mult,mult)
variable niu(n+m,n)
variable W(m+n,n)

expression H2n

% Define Phix, Phiu and Psi
for ind = 1:mult
    Phix{ind} = phix(:,:,ind);
    Phiu{ind} = phiu(:,:,ind);
    Psi{ind} = [Q^0.5*phix(:,:,ind);R^0.5*phiu(:,:,ind)];
end


for ind = mult:-1:1
    if ind<mult
        C_stack(:,:,mult+1-ind) = Psi{ind}-Psi{ind+1};
    else
        C_stack(:,:,mult+1-ind) = [Psi{ind}];
    end
end

for ind1 = 1:m+n
    for ind2 = 1:n
        Cr{ind1,ind2} = squeeze(C_stack(ind1,ind2,:))';
        Pr{ind1,ind2} = squeeze(P(ind1,ind2,:,:));
    end
end

H2n = ones(1,n+m)*niu*ones(n,1);

minimize H2n
subject to

% Compute H2 norm, we compute it element wise
for ind1 = 1:m+n
    for ind2 = 1:n
        tP = Pr{ind1,ind2};
        tC = Cr{ind1,ind2};
        tP > 0;
        [Ar'*tP+tP*Ar,tP*Br; Br'*tP,-1]<=0
        [tP,tC';tC,W(ind1,ind2)]>=0
        W(ind1,ind2)<=niu(ind1,ind2)
    end
end


% Feasibility constraints
Phix{1} == eye(n);
for ind = 1:mult-1
   norm(Phix{ind+1}+(p*eye(n)-A)*Phix{ind}-B*Phiu{ind}) <= 1e-7
end
norm((p*eye(n)-A)*Phix{mult}-B*Phiu{mult}) <= 1e-7


cvx_end


%% Phix and Phiu in TF, verify if the feasibility constraints hold
Phix_tf = zeros(n);
Phiu_tf = zeros(m,n);
Psi_tf = zeros(n+m,n);
for ind = 1:mult
    Phix_tf = Phix_tf+Phix{ind}/(s-p)^ind;
    Phiu_tf = Phiu_tf+Phiu{ind}/(s-p)^ind;
    Psi_tf = Psi_tf+[Q^0.5*Phix{ind};R^0.5*Phiu{ind}]/(s-p)^ind;
end

err = (s*eye(n)-A)*Phix_tf-B*Phiu_tf-eye(n);
% e_norm = norm(err);

% Forming A,B,C for norm computing
Am = zeros(n*mult);
Bm = zeros(n*mult,n);

for ind = 1:mult
    Bm((ind-1)*n+1:ind*n,:) = eye(n);
    if ind<mult
        Am(1+(ind-1)*n:ind*n,1+(ind-1)*n:(1+ind)*n) = [p*eye(n),eye(n)];
    else
        Am(1+(ind-1)*n:ind*n,1+(ind-1)*n:(ind)*n) = p*eye(n);
    end
end
Cm = [];
for ind = 1:mult
    Cm = [Cm,C_stack(:,:,ind)];
end

%First H2 norm verification:
clsys = ss(Am,Bm,Cm,0);
% H2cl = norm(clsys,2)
% H2cl2 = norm(Psi_tf,2)
H2cvx = H2n^0.5;
H2_normalized = H2cvx/(H2norm_LQR)

Phix;
Phiu;

%% Save data

