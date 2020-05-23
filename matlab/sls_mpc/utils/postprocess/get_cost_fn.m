function obj = get_cost_fn(params, x, u)

obj=0;
for t=1:params.tHorizon_
    obj = obj + x(:,t)'*params.Q_*x(:,t)+u(:,t)'*params.R_*u(:,t);
end
obj = obj + x(:,t+1)'*params.Q_*x(:,t+1);

end