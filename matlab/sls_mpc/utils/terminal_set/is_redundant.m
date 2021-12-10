function isRedundant = is_redundant(H, h, i)
% Check whether row i of Hx <= h  is redundant

h_alt    = h;
h_alt(i) = h(i) + 1;

% Note: there is a more efficient way via Clarkson/RayShoot but it is
% a little more involved to code so we are doing it the naive way

% Set up LP
% z = max H(i,:)*x
%     s.t. H(not_i, :)*x <= h(not_i, :)
%          H(i)*x        <= h(i) + 1
% z > h(i) indicates redundancy

model.obj        = -H(i,:);
model.A          = sparse(H);
model.rhs        = h_alt;
model.lb         = -inf*ones(size(H,2),1); % default is 0

gParams.outputflag     = 0;
gParams.dualReductions = 0; % clearly output unbounded or not
gParams.numericFocus   = 3; % avoid numerical issues
result = gurobi(model, gParams);

% Adjust EPS based on A size; more tolerance for very large A
EPS = 1e-8;
if size(model.A, 1) > 40
    EPS = 1e-5;
elseif size(model.A, 1) > 30
    EPS = 1e-6;
elseif size(model.A, 1) > 20
    EPS = 1e-7;
end

isRedundant = 0;
if strcmp(result.status, 'NUMERIC')
    mpc_error('Gurobi reports numeric issues; model may be too large')
else
    isRedundant = (strcmp(result.status, 'UNBOUNDED')) || (-result.objval <= h(i) + EPS);
end

if ~isRedundant
    size(model.A);
end
end
