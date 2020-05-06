function fzuw = f_glyco_log(z, u, w, q, k, a)

z1=z(1); z2=z(2);
u1=u(1); u2=u(2);

t1 = 2*exp(a*z2-z1+u1);
t2 = 2*k*exp(u2);
t3 = 2*q*exp((a-1)*z2+u1);
t4 = 2*k*(q+1)*exp(z1-z2+u2);

fzu = [ t1 - t2;
       -t3 + t4];
   
e2  = [0; exp(-z2)];

fzuw = fzu + e2*w;

end