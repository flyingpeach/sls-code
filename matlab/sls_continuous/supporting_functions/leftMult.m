function prod = leftMult(Phiu,B)
n = size(Phiu{1},2);
Polenum = size(Phiu,1);

cvx_begin quiet
expression prod(n,n*Polenum);
 
for l = 1:Polenum
    prod(:,(l-1)*n+1:l*n) = B*Phiu{l};
end

cvx_end
end