function fxuw = f_bilo(x, u, q, w)

x1=x(1); x2=x(2); x4=x(4);
u1=u(1); u2=u(2); u3=u(3);

v1 = x1*x4*exp(u1);
v2 = x2*x4*exp(u2);
v3 = x4*exp(u3);

s1 = [-1 1 0 -q 0]';
s2 = [0 -1 1 (q+1) 0]';
s3 = [0 0 0 -1 1]';

v = [v1; v2; v3];
S = [s1 s2 s3];

fxuw = S*v + w;

end