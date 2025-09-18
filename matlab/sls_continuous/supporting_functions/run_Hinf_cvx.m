function [Hinfnorm,Psi,Phix,Phiu,phix_supp,phiu_supp,status] = run_Hinf_cvx(sysA,sysB,Q,R,Poles,flag,d)
[n,m] = size(sysB);
K = size(Poles,1);

% State space realization of Phi
A = zeros(n*K);
B = zeros(n*K,n);
for j = 1:K
    A(j*n-n+1:j*n,j*n-n+1:j*n) = Poles(j)*eye(n);
    B(j*n-n+1:j*n,:) = eye(n);
end
% Creating a real realization
T = [];
U = [(1+1i)/2*eye(n),(1-1i)/2*eye(n);(1-1i)/2*eye(n),(1+1i)/2*eye(n)];
index=1;
while (index<=K)
    if(imag(Poles(index))~=0)
        T = blkdiag(T,U);
        index = index+2;
    else
        T = blkdiag(T,eye(n));
        index = index+1;
    end
end
A0 = T'*A*T;
B0 = T'*B;

cvx_begin sdp quiet%%

%Expressions for phix and phiu
variable X(n*K,n*K) hermitian
variable Gamma

if flag==0
    variable tPhix(n,n*K)
    variable tPhiu(m,n*K)
    phix_supp = ones(n);
    phiu_supp = ones(m,n);
else% Sparse Hinf
    PhixSupp=zeros(n,n*K);
    PhiuSupp=zeros(m,n*K);
    %Computing the support of the variables
    Bsupp = sysB~=0;
    if flag==1
        Acomm = adjust_grid_sys_locality(sysA);
        Asupp = Acomm~=0;
    else
        Asupp = sysA~=0;
    end
    phix_supp = (Asupp+eye(n))^d ~=0;
    phiu_supp = Bsupp'*phix_supp ~=0;
    for index = 1:K
        PhixSupp(:,(index-1)*n+1:index*n) = phix_supp;
        PhiuSupp(:,(index-1)*n+1:index*n) = phiu_supp;
    end
    Tsupp = T~=0;
    tPhix_supp = PhixSupp*Tsupp ~= 0;
    tPhiu_supp = PhiuSupp*Tsupp ~= 0;
    dimPhix = sum(tPhix_supp(:));
    dimPhiu = sum(tPhiu_supp(:));

    variable vPhix(dimPhix,1)
    variable vPhiu(dimPhiu,1)
    expression tPhix(n,n*K)
    expression tPhiu(m,n*K)
    tPhix(tPhix_supp) = vPhix;
    tPhiu(tPhiu_supp) = vPhiu;
end
C0 = [Q^(1/2)*tPhix;R^(1/2)*tPhiu];

%set system C matrix
J1 = [A0'*X'+X*A0,X*B0,C0';B0'*X',-Gamma*eye(n),zeros(n,n+m);
        C0,zeros(n+m,n),-Gamma*eye(n+m)];

minimize(Gamma)
subject to 
-J1 == semidefinite(n*(K+2)+m)
% tPhix*A0-sysA*Q^(1/2)*tPhix == sysB*tPhiu %This is (p_l*I-sysA)*Phix == sysB*Phiu
norm(tPhix*B0-eye(n)) <= 1e-14
norm(tPhix*A0-sysA*Q^(1/2)*tPhix-sysB*tPhiu)  <= 1e-14
% tPhix*B0-eye(n)==0
% tPhix*A0-sysA*Q^(1/2)*tPhix-sysB*tPhiu==0

cvx_end
status = cvx_status;

%The H_infty norm
Hinfnorm = Gamma;

Phix = cell(K,1);
Phiu = cell(K,1);
Psi = cell(K,1);

for index = 1:K
    Phix{index} = tPhix(:,(index-1)*n+1:index*n);
    Phiu{index} = tPhiu(:,(index-1)*n+1:index*n);
    Psi{index} = [Q^0.5*Phix{index};R^0.5*Phiu{index}];
end

end