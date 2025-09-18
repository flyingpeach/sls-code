%% This is the function that computes the cost 
function [H2norm,Psi] = computeH2norm(Phix,Phiu,Q,R,Poles,pij)
% Get parameters related to dimension
% pu1 = Phiu{1};
[m,n] = size(Phiu{1});
Polenum = size(Phiu,1);

% Compute the matrix M used for cost computation of each element
M = zeros(Polenum);
for p = 1:Polenum
    for j = 1:Polenum
        M(p,j) = -1/(conj(Poles(p))+Poles(j));
    end
end
sM = M^(0.5);

Psi = cell(Polenum,1);

sQ = Q^0.5;
sR = R^0.5;
for l = 1:Polenum
    Psi{l} = [sQ*Phix{l};sR*Phiu{l}];
end

H2norm = 0;

for p = 1:m+n
    for j = 1:n
        for l = 1: Polenum
            Pl = Psi{l};
            pij(l) = Pl(p,j);
        end
        vect = sM*pij;
        H2norm = H2norm + (vect')*vect;
    end
end
end