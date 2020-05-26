function plot_b(x_series, times, timeCents, iters, localityBool)

if localityBool
    xStr = '$$\#\ subsystems\ in\ localized\ region$$';
else
    xStr = '$$\#\ pendulums\ in\ network$$';
end

figure()
subplot(1,2,1)
hold on
plot(x_series, times,'m-s','LineWidth',2)
plot(x_series, timeCents,'b-s','LineWidth',2)  
xlabel(xStr, 'Interpreter','latex','Fontsize', 10)
ylabel('$$Avg\ runtime\ per\ MPC\ iteration\ for\ each\ state\ (s)$$','Interpreter','latex','Fontsize', 10)
leg1 = legend('$$Localized\ MPC\ using\ ADMM$$', '$$Centralized\ MPC$$');
set(leg1,'Interpreter','latex','Fontsize', 8);

subplot(1,2,2)
plot(x_series, iters,'m-s','LineWidth',2)
xlabel(xStr, 'Interpreter','latex','Fontsize', 10)
ylabel('$$Avg\ \#\ ADMM\ iters\ per\ MPC\ iteration\ for\ each\ state\ (s)$$','Interpreter','latex','Fontsize', 10)

end