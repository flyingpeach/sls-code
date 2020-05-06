function [A, B] = linearize_glyco_log(z_, u_, q, k, a)
% get matrices for linearized dynamics
% linearize about z_, u_; don't forget to add f(z_, u_, w_)
 
z1=z_(1); z2=z_(2);
u1=u_(1); u2=u_(2);

t1 = 2*exp(a*z2-z1+u1);
t2 = 2*k*exp(u2);
t3 = 2*q*exp((a-1)*z2+u1);
t4 = 2*k*(q+1)*exp(z1-z2+u2);

A = [-t1             a*t1;
      t4   -(a-1)*t3 - t4];

B = [ t1 -t2;
     -t3  t4];

end