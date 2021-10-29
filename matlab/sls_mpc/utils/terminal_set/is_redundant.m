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

gParams.outputflag = 0;
result = gurobi(model, gParams);

EPS = 1e-8;
isRedundant = -result.objval <= h(i) + EPS;

end
