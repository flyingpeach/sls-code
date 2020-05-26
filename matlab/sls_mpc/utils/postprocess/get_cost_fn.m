function obj = get_cost_fn(params, x, u)

obj = 0;

tHorizon = size(x, 2);

Q = params.QSqrt_*params.QSqrt_;
R = params.RSqrt_*params.RSqrt_;

for t=1:tHorizon-1
    obj = obj + x(:,t)'*Q*x(:,t)+u(:,t)'*R*u(:,t);
end

t   = tHorizon; % the last u is not calculated, but the last x is 
obj = obj + x(:, t)'*Q*x(:, t);

end