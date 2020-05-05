function [A, B] = discretize(Ac, Bc, Ts)

Nx = size(Ac, 1);
A  = eye(Nx) + Ac*Ts;
B  = Bc*Ts;

end