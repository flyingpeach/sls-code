function [M1, M2, MSum, MbSum] = coupled_row_setup(x_loc, y, z, selfIdx)

n   = length(x_loc);
nc  = length(z); % # of coupled subsystems (including self)

MSum  = 0; % "Quadratic" terms for the mu/2 objective term
MbSum = 0; % "Linear" terms for the mu/2 objective term
for j = 1:nc
    if j < selfIdx
            Mj      = zeros(1,nc+n-1); 
            Mj(n+j) = 1;
    elseif j == selfIdx
            Mj      = [x_loc' zeros(1,nc-1)]; % Implicitly enforce Xi = Phi*x0
    elseif j > selfIdx
            Mj        = zeros(1,nc+n-1); 
            Mj(n+j-1) = 1;
    end

    MSum  = MSum + (Mj'*Mj);
    b     = z{j} - y{j};    
    MbSum = MbSum + Mj'*b;
end

% These terms are not related to the mu/2 objective terms
myEye = eye(nc-1);
M1    = [eye(n) zeros(n,nc-1)];
M2    = [zeros(selfIdx-1,n)  myEye(1:selfIdx-1,:); 
         x_loc'            zeros(1,nc-1); 
         zeros(nc-selfIdx,n) myEye(selfIdx:end,:)];                                 
                                 
end