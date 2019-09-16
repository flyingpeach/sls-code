%% sanity checks: R, M should be solutions for Rc, Mc for Tc >= T
setup2;

% for rank/zero conditions, try to match the precision of cvx_precision low
% http://cvxr.com/cvx/doc/solver.html#solver-precision
eps = 2.22e-16;
tol = eps.^(3/8);

close all;

Rs = slsOuts.R_;
Ms = slsOuts.M_;

T = slsParams.tFIR_;
Tc = T;

%% test: does R, M solve equation? (block)

% zero padding for help
for t=T+1:Tc
    Rs{t} = zeros(sys.Nx, sys.Nx);
    Ms{t} = zeros(sys.Nu, sys.Nx);
end

normBlock = 0;

Dellc{1}    = eye(sys.Nx);
for t=2:Tc
    Dellc{t} = Rs{t} - sys.A*Rs{t-1} - sys.B2*Ms{t-1};
end
Dellc{Tc+1} = -sys.A*Rs{Tc} - sys.B2*Ms{Tc};
for t=Tc+2:Tc+T
    Dellc{t} = zeros(sys.Nx, sys.Nx); % zero padding for ease of calculation
end

% enforce LHS = RHS constraints
for t=1:Tc
    Rsum = zeros(sys.Nx, sys.Nx); Msum = zeros(sys.Nu, sys.Nx);
    for k=1:min(T, t) % convolve
        Rsum = Rsum + Rs{k} * Dellc{t-k+1};
        Msum = Msum + Ms{k} * Dellc{t-k+1};
    end
    normBlock = normBlock + norm(Rs{t} - Rsum); % Rs{t} == Rsum
    normBlock = normBlock + norm(Ms{t} - Msum);
end
for t=Tc+1:Tc+T
    Rsum = zeros(sys.Nx, sys.Nx); Msum = zeros(sys.Nu, sys.Nx);
    for k=1:min(T, t)
         Rsum = Rsum + Rs{k} * Dellc{t-k+1};
         Msum = Msum + Ms{k} * Dellc{t-k+1};
    end
    normBlock = normBlock + norm(Rsum); % Rsum == 0
    normBlock = normBlock + norm(Msum);
end

normBlock

%% test: does R, M solve equation? (precise)
Rs_ = zeros(sys.Nx * (Tc-1), sys.Nx); % Rc{1} == 1 so it's not a variable
Ms_ = zeros(sys.Nu * Tc, sys.Nx);

% helper function to convert cells of block matrices into giant matrix
% copied from get_F
get_range = @(idx, size) (size*(idx-1)+1:size*(idx-1)+size);

for t = 1:T
    tx    = get_range(t-1, sys.Nx);
    tu    = get_range(t, sys.Nu);
    if t > 1
        Rs_(tx,:) = Rs{t};
    end
    Ms_(tu,:) = Ms{t};
end

% constraints
F  = get_F(sys, slsParams, slsOuts, Tc);
F1 = F(:, 1:sys.Nx);
F2 = F(:,sys.Nx+1:end);

expectZeroPrecise = F2 * [Rs_; Ms_] + F1;
normPrecise       = norm(expectZeroPrecise)

