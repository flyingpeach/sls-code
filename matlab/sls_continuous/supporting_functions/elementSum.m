function phix_sum = elementSum(Phix)
K = size(Phix,1);
phix_sum = Phix{1};
for p = 2:K
    phix_sum = phix_sum+Phix{p};
end
end