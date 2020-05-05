function fxuw = f_bilo(x, u, w, q, k)

x1=x(1); x2=x(2); x4=x(4);
u1=u(1); u2=u(2); u3=u(3);
k1=k(1); k2=k(2); k3=k(3);

S = [-1    0   0;
      1   -1   0; 
      0    1   0;
     -q  q+1  -1;
      0    0   1];

v = [k1*x1*x4*exp(u1);
     k2*x2*x4*exp(u2);
     k3*x4*exp(u3)];

e1   = [1 0 0 0 0]';
fxuw = S*v + e1*w;

end