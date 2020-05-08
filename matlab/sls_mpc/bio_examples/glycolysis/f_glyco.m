function fxuw = f_glyco(x, u, w, q, k, a)

Q = [ 1   -1  0;
     -q  q+1 -1
      0    0  1];

K = diag(k);

v = [x(2).^a * exp(u(1));
     x(1)    * exp(u(2));
     x(2)    * exp(u(3))];
 
e2 = [0 1 0]';

fxuw = Q*K*v + e2*w;

end