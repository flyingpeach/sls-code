function phi_x = h2_sls(sys, params)
% INFINITE HORIZON LOCALIZED LQR
% Author: Jing Yu (&edited by Lisa Li)
%
% This code performs the computations in eq.(14) from the article
% "Localized and Distributed H2 State Feedback Control" by J. Yu, Y.-S. Wang
% and J. Anderson. ArXiv: https://arxiv.org/abs/2010.02440
%
% sys    : LTISystem object; we will use A and B2
% params : MPCParams object; we will use locality_, QSqrt_, and RSqrt_
eps = 1e-8;

%% Setup and allocate variables
Nx = sys.Nx; A = sys.A; B = sys.B2;
d  = params.locality_;

% Placeholder to transform R to Q-size
R2Q = eye(Nx);

xSparsity       = abs(A^(d-1)) > 0;
xSparsityExt    = abs(A^(d)) > 0;
boundaryPattern = abs(xSparsityExt - xSparsity) > 0;
uSparsity       = (abs(B') > 0) * abs(A^(d)) > 0;

K  = cell(Nx, 1);
S  = cell(Nx, 1);

% These are not the standard closed-loop maps; rather they are
% the multiplier terms used to recursively calculate the infinite-horizon
% closed-loop maps as per eq (13) in https://arxiv.org/pdf/2010.02440
phi_x = cell(Nx, 1);
phi_u = cell(Nx, 1);

%% Local K Computation
for i = 1:Nx
    internal = find(xSparsity(:,i)); % find which states are internal
    boundary = find(boundaryPattern(:,i)); % find which states are on the boundary
    control  = find(uSparsity(:,i));
        
    A_nn = A(internal, internal);
    A_bn = A(boundary, internal);
    
    B_n  = B(internal, control);
    B_b  = B(boundary, control);
    
    R2Q_ = R2Q(internal, control);
     
    lqr_A = A_nn - B_n*pinv(full(B_b))*A_bn;
    lqr_B = B_n*(eye(size(B_n,2)) - pinv(full(B_b))*B_b);

    Q_loc = params.QSqrt_(xSparsity(:, i), xSparsity(:, i));
    R_loc = params.RSqrt_(uSparsity(:, i), uSparsity(:, i));
    
    Q_ = Q_loc - R2Q_*R_loc*pinv(full(B_b))*A_bn;
    R_ = R_loc * (eye(size(B_b,2)) - pinv(full(B_b))*B_b);
    
    lqr_Q = Q_'*Q_;
    lqr_R = R_'*R_;
    
    % Slight hack: add nonzero element to ensure positive definite
    lqr_R = lqr_R + eps*eye(size(lqr_R));
    
    [K{i}, S{i}, ~] = dlqr(lqr_A, lqr_B, lqr_Q, lqr_R);
    phi_u{i} = -pinv(full(B_b))*A_bn - (eye(size(B_n,2)) - pinv(full(B_b))*B_b)*K{i};
    phi_x{i} = lqr_A - lqr_B*K{i};
    
end

end
