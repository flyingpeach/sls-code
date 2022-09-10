function [Phi_x, Phi_u, nomCost] = phi_step(sys, params, D, beta)
% TODO: this can be combined with existing sls_base code by taking
% advantage of the addConstraint and addObjective methods

A = sys.A; B = sys.B2; 
T = params.tFIR_;

QSqrt = params.QSqrt_; RSqrt = params.RSqrt_;
Nx = size(B,1); Nu = size(B,2); 

% Translate A into adjacency matrix, with self-connections
adj = ((A ~= 0) + eye(Nx) + (A' ~= 0)) ~= 0;

% Support sizes (based solely on local patch size d)
PhixSupp = adj^(params.locality_-1) ~= 0;
PhiuSupp = ((B ~=0)' * PhixSupp) ~= 0;

suppSizePhix = sum(sum(PhixSupp));
suppSizePhiu = sum(sum(PhiuSupp));
numSupps     = suppSizePhix*T + suppSizePhiu*T; % number nonzero elements

cvx_begin quiet
expression Phi_xs(Nx, Nx, T) % 3D matrices
expression Phi_us(Nu, Nx, T)
variable PhiSupps(numSupps)

% Populate decision variables for ease-of-use
Phi_x = cell(T,1); Phi_u = cell(T,1); % List of 2D matrices
for k=1:T
    Phi_x{k} = Phi_xs(:,:,k); Phi_u{k} = Phi_us(:,:,k);
end

spot = 0;
for k=1:T % Populate expressions with support values
    Phi_x{k}(PhixSupp) = PhiSupps(spot+1:spot+suppSizePhix);
    spot = spot + suppSizePhix;
    
    Phi_u{k}(PhiuSupp) = PhiSupps(spot+1:spot+suppSizePhiu);
    spot = spot + suppSizePhiu;
end

% SLS feasibility constraints
Phi_x{1} == eye(Nx);
Phi_x{T} == zeros(Nx, Nx);
for k=1:T-1
    Phi_x{k+1} == A*Phi_x{k} + B*Phi_u{k};
end

nomCost = 0;
for k=1:T
    % Nominal performance objective (LQR)
    vectNom = vec(blkdiag(QSqrt, RSqrt) * [Phi_x{k}; Phi_u{k}]);
    nomCost = nomCost + vectNom'*vectNom;    
end

M        = make_m(params, Phi_x, Phi_u);
robBound = get_bound(D, M, params.stabType_);

% Robust stability progress
if ~isinf(beta)
    robBound <= beta;
end

minimize(nomCost);
cvx_end

end