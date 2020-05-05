function [A, B] = linearize_glyco(x_, u_, q, k, a)
% get matrices for linearized dynamics
% linearize about x_, u_; don't forget to add f(x_, u_, w_)
 
x1 = x_(1); x2 = x_(2);
u1 = u_(1); u2 = u_(2);

eu1 = exp(u1); eu2 = exp(u2);

Q = [ 1   -1;
     -q  q+1];

qa = [      0  2*a*x2.^(a-1)*eu1;
      2*k*eu2                 0];

qb = [2*x2.^a*eu1          0;
                0 2*k*x1*eu2];

A = Q*qa;
B = Q*qb;

end