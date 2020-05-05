function [A, B] = linearize_bilo(x_, u_, q, k)
% get matrices for linearized dynamics
% linearize about x_, u_; don't forget to add f(x_, u_, w_)

S = [-1    0   0;
      1   -1   0; 
      0    1   0;
     -q  q+1  -1;
      0    0   1];
  
x1 = x_(1); x2 = x_(2); x4 = x_(4);
u1 = u_(1); u2 = u_(2); u3 = u_(3);
k1 = k(1) ; k2 = k(2) ; k3 = k(3);

eu1 = exp(u1); eu2 = exp(u2); eu3 = exp(u3);


sa = [k1*x4*eu1          0  0  k1*x1*eu1  0;
              0  k2*x4*eu2  0  k2*x2*eu2  0;
              0          0  0     k3*eu3  0];
  
sb = [k1*x1*x4*eu1             0          0;
                 0  k2*x2*x4*eu2          0;       
                 0             0  k3*x4*eu3];

A = S*sa;
B = S*sb;

end