function obj = get_cost_fn(Q, S, tSim, x, u)

obj=0;
for t=1:tSim
    obj = obj + x(:,t)'*Q*x(:,t)+u(:,t)'*S*u(:,t);
end
obj = obj + x(:,t+1)'*Q*x(:,t+1);

end