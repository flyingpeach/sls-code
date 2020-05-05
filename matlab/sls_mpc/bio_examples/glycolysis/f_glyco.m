function fxuw = f_glyco(x, u, w, q, k, a)

Q = [ 1   -1;
     -q  q+1];

v = [2*x(2).^a*exp(u(1));
     2*k*x(1)*exp(u(2))];

e2 = [0; 1];

fxuw = Q*v + e2*w;

end