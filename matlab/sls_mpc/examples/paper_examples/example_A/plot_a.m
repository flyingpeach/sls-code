function plot_a(params, x, u, xVal, uVal, myTitle)

% Calculate costs + plot 
time   = 1:size(x,2);
obj    = get_cost_fn(params, x, u);
objVal = get_cost_fn(params, xVal, uVal);

% Print costs (sanity check: should be close)
fprintf('Distributed cost: %f\n', obj);
fprintf('Centralized cost: %f\n', objVal);

figure()
plot(time,xVal(1,:),'b',time,x(1,:),'*b',time,xVal(3,:),'g',time,x(3,:),'*g')
xlabel('$$Time$$','interpreter','latex','Fontsize', 10)
ylabel('$$\theta_{1},\ \theta_{2}$$','Interpreter','Latex','Fontsize', 10)
leg1 = legend('$$\theta_{1}\ Centralized\ MPC$$', '$$\theta_{1}\ Localized\ MPC\ using\ ADMM$$','$$\theta_{2}\ Centralized\ MPC$$', '$$\theta_{2}\ Localized\ MPC\ using\ ADMM$$');
set(leg1,'Interpreter','latex'); set(leg1, 'Fontsize', 8)
title(strcat(myTitle, ', Subsystems 1 and 2'));

end
