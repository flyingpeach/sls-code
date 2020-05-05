function [A, B] = linearize_bilo(x_, u_, q, k)
% get matrices for linearized dynamics
% linearize about x_, u_; don't forget to add f(x_,u_, w)

x1 = x_(1); x2 = x_(2); x4 = x_(4);
u1 = u_(1); u2 = u_(2); u3 = u_(3);
k1 = k(1) ; k2 = k(2) ; k3 = k(3);

eu1 = exp(u1); eu2 = exp(u2); eu3 = exp(u3);

a1 = [-k1*x4*eu1, 0, 0, -k1*x1*eu1, 0];
a2 = [k1*x4*eu1, -k2*x4*eu2, 0, (k1*x1*eu1 - k2*x2*eu2), 0];
a3 = [0, k3*x4*eu2, 0, k2*x2*eu2, 0];
a4 = [-q*k2*x4*eu1, (q+1)*k2*x4*eu2, 0, (-q*k1*x1*eu1 + (q+1)*k2*x2*eu2 - k3*eu3), 0];
a5 = [0, 0, 0, k3*eu3, 0];

b1 = [-k1*x1*x4*eu1, 0, 0];
b2 = [k1*x1*x4*eu1, -k2*x2*x4*eu2, 0];
b3 = [0, k2*x2*x4*eu2, 0];
b4 = [-q*k1*x1*x4*eu1, (q+1)*k2*x2*x4*eu2, -k3*x4*eu3];
b5 = [0, 0, k3*x4*eu3];

A = [a1; a2; a3; a4; a5];
B = [b1; b2; b3; b4; b5];

end