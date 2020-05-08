function [A, B] = linearize_glyco(x_, u_, q, k, a)
% get matrices for linearized dynamics
% linearize about x_, u_; don't forget to add f(x_, u_, w_)
 
x1 = x_(1); x2 = x_(2); x3 = x_(3);
u1 = u_(1); u2 = u_(2); u3 = u_(3);

eu1 = exp(u1); eu2 = exp(u2); eu3 = exp(u3);

Q = [ 1   -1  0;
     -q  q+1 -1
      0    0  1];

K = diag(k);

qa = [  0  a*x2.^(a-1)*eu1  0;
      eu2                0  0
        0              eu3  0];

qb = [x2.^a*eu1       0       0;
              0  x1*eu2       0
              0       0  x2*eu3];

A = Q*K*qa;
B = Q*K*qb;

end