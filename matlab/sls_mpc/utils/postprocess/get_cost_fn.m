function obj = get_cost_fn(params, x, u)

obj=0;

Q = params.QSqrt_*params.QSqrt_;
R = params.RSqrt_*params.RSqrt_;

for t=1:params.tHorizon_
    obj = obj + x(:,t)'*Q*x(:,t)+u(:,t)'*R*u(:,t);
end
obj = obj + x(:,t+1)'*Q*x(:,t+1);

end