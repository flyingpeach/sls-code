function [H2norm,Psi,Phix,Phiu,phix_supp,phiu_supp,status] = run_H2_cvx(A,B,Q,R,Poles,flag,d)

warning('off');
% if(sparsity == 'grid')
% flag = 1;
% elseif(sparsity == 'chain')
% flag = 2;
% else % 'dense' is considered 
% flag = 0;
% end

n = size(A,1); %dimension of the system
m = size(B,2);
Polenum = size(Poles,1)

cvx_begin quiet

expression Psi_in(n+m,n);
expression pij(Polenum,1);

%These are for non-sparse cvx
if(flag==0) %Dense cvx
    variable phi_x(n,n,Polenum) complex;
    variable phi_u(m,n,Polenum) complex;
    Phix = cell(Polenum,1);
    Phiu = cell(Polenum,1);
    for index = 1:Polenum
        Phix{index} = phi_x(:,:,index);
        Phiu{index} = phi_u(:,:,index);
    end
    phix_supp = ones(n);
    phiu_supp = ones(m,n);

else
    expression phi_x(n,n,Polenum);
    expression phi_u(m,n,Polenum);
    %Computing the support of A, B, Phix and Phiu
    Bsupp = B~=0;

    if flag==1 
        % flag = 1, Grid Model
        Acomm = adjust_grid_sys_locality(A);
        Asupp = Acomm~=0;
    else
        % flag = 2, Chain Model
        Asupp = A~=0;
    end
    phix_supp = (Asupp+eye(n))^d ~=0;
    phiu_supp = Bsupp'*phix_supp ~=0;
    Phix = cell(Polenum,1);
    Phiu = cell(Polenum,1);
    
    %Create Cell variables for Phix and Phiu
    for index = 1:Polenum
        Phix{index} = phi_x(:,:,index);
        Phiu{index} = phi_u(:,:,index);
    end
    dim_phix_vec = sum(phix_supp(:));
    dim_phiu_vec = sum(phiu_supp(:));
    variable vphi_x(dim_phix_vec,Polenum) complex;
    variable vphi_u(dim_phiu_vec,Polenum) complex;
    for index = 1:Polenum
        Phix{index}(phix_supp) = vphi_x(:,index);
        Phiu{index}(phiu_supp) = vphi_u(:,index);
    end
end

[H2norm,Psi] = computeH2norm(Phix,Phiu,Q,R,Poles,pij);
s_phi_x = elementSum(Phix);
p_phi_x = pI_m_A(Phix,Poles,A);
b_phi_u = leftMult(Phiu,B);

minimize(H2norm);

subject to
% 
s_phi_x-eye(n)==0
p_phi_x-b_phi_u==0

% norm(s_phi_x- eye(n),2) <= 1e-6;
% s_phi_x- eye(n) <= (1e-5)*ones(n);
% 
% s_phi_x- eye(n) >= -(1e-5)*ones(n);
% norm(p_phi_x-b_phi_u,2) <= 1e-6;

index=1;
while index<= Polenum
    if imag(Poles(index)) == 0
        % norm(imag(phi_x(:,:,index)),2) <= 1e-6;
        imag(phi_x(:,:,index)) <= (1e-6)*ones(n);
        imag(phi_x(:,:,index)) >= -(1e-6)*ones(n);
        % index = index+1;
    else
        % norm(real(phi_x(:,:,index))-real(phi_x(:,:,index+1)),2) <= 1e-6;
        % norm(imag(phi_x(:,:,index))+imag(phi_x(:,:,index+1)),2) <= 1e-6;
        real(phi_x(:,:,index))-real(phi_x(:,:,index+1)) <= (1e-10)*ones(n);
        real(phi_x(:,:,index))-real(phi_x(:,:,index+1)) >= -(1e-10)*ones(n);
        imag(phi_x(:,:,index))+imag(phi_x(:,:,index+1)) <= (1e-10)*ones(n);
        imag(phi_x(:,:,index))+imag(phi_x(:,:,index+1)) >= -(1e-10)*ones(n);
        index = index+2;
    end
end

cvx_end

status = cvx_status;

% index=1;
% % Removing some annoying terms
% while index<=Polenum
%     if imag(Poles(index)) == 0
%         Psi{index} = real(Psi{index});
%         phi_x(:,:,index) = real(phi_x(:,:,index));
%         phi_u(:,:,index) = real(phi_u(:,:,index));
%         index = index+1;
%     else
%         Psi{index+1} = conj(Psi{index});
%         phi_x(:,:,index+1) = conj(phi_x(:,:,index));
%         phi_u(:,:,index+1) = conj(phi_u(:,:,index));
%         index = index+2;
%     end
% end
% 
% %Trying out another way to compute the norms
% modelA = zeros(n*Polenum);
% modelB = zeros(n*Polenum,n);
% modelC = zeros(m+n,n*Polenum);
% for index = 1:Polenum
%     modelA((index-1)*n+1:index*n,(index-1)*n+1:index*n) = Poles(index)*eye(n);
%     modelB((index-1)*n+1:index*n,:) = eye(n);
%     modelC(:,(index-1)*n+1:index*n) = Psi{index};
% end
% 
% ss_model = ss(modelA,modelB,modelC,0);
% 
% H2_ss = norm(ss_model,2)
% H2norm
end