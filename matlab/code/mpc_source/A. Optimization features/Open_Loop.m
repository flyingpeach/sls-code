setup_system_a;

x_VAL(:,1) = x0;
for t = 1:tSim
    x_VAL(:,t+1) = A*x_VAL(:,t);
end

figure(2)
plot(1:tSim+1,x_VAL(1,:),1:tSim+1,x_VAL(3,:))