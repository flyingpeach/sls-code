function [flag,inds] = Hinf_sparsity_sanity_check(Phix,Phiu,phix_supp,phiu_supp)
inds = [];
flag = 1;

Polenum = size(Phix,1);

for index = 1:Polenum
    nPhix = Phix{index};
    nPhiu = Phiu{index};
    bPhix = (nPhix~=0)&(~phix_supp);
    bPhiu = (nPhiu~=0)&(~phiu_supp);
    temp = sum(bPhiu(:))+sum(bPhix(:));
    if temp ~=0
        inds = [inds;index];
        flag = 0;
    end
end

end