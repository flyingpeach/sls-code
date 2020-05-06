h = 3;
g = 0.1;

cpfk      = xs(2,:).^a ./ (1 + xs(2,:).^(2*h));
cpfk_mine = xs(2,:).^a .* exp(us(1,:));

cpk       = xs(1,:) ./ (1 + xs(2,:).^(2*g));
cpk_mine  = xs(1,:) .* exp(us(2,:));

figure();
subplot(2,1,1);
hold on;
plot(cpfk)
plot(cpfk_mine)
legend('ref','mine');

subplot(2,1,2);
hold on;
plot(cpk)
plot(cpk_mine)