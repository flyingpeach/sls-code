function [M1, M2, MSum, MbSum] = coupled_row_setup(x_loc, y, z, sIdx)

n   = length(x_loc);
nc  = length(z); % # of coupled subsystems (including self)

MSum  = 0; % "Quadratic" terms for the mu/2 objective term
MbSum = 0; % "Linear" terms for the mu/2 objective term
b     = z - y;
for j = 1:nc
    if j < sIdx
            Mj      = zeros(1,nc+n-1); 
            Mj(n+j) = 1;
    elseif j == sIdx
            Mj      = [x_loc' zeros(1,nc-1)]; % Implicitly enforce Xi = Phi*x0
    elseif j > sIdx
            Mj        = zeros(1,nc+n-1); 
            Mj(n+j-1) = 1;
    end

    MSum  = MSum + (Mj'*Mj);
    MbSum = MbSum + Mj'*b(j);
end

% These terms are not related to the mu/2 objective terms
myEye = eye(nc-1);
M1    = [eye(n) zeros(n,nc-1)];
M2    = [zeros(sIdx-1,n)  myEye(1:sIdx-1,:); 
         x_loc'            zeros(1,nc-1); 
         zeros(nc-sIdx,n) myEye(sIdx:end,:)];                                 
                                 
end