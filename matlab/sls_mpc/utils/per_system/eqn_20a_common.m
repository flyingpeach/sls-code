function [M1, M2, MSum, MbSum] = eqn_20a_common(x_loc, y_, z_, selfIdx)

n   = length(x_loc);
nc  = length(z_);    % # of coupled subsystems (including self)

MSum  = 0; % "Quadratic" terms for the mu/2 objective term
MbSum = 0; % "Linear" terms for the mu/2 objective term
for j = 1:nc
    if j < selfIdx
            Mj      = zeros(1,nc+n-1); 
            Mj(n+j) = 1;
    elseif j == selfIdx
            Mj      = [x_loc' zeros(1,nc-1)];
    elseif j > selfIdx
            Mj        = zeros(1,nc+n-1); 
            Mj(n+j-1) = 1;
    end

    MSum  = MSum + (Mj'*Mj);

    b     = z_{j} - y_{j};    
    MbSum = MbSum + Mj'*b;
end

% These terms are not related to the mu/2 objective terms
myEye = eye(nc-1);
M1    = [eye(n) zeros(n,nc-1)];
M2    = [zeros(selfIdx-1,n)  myEye(1:selfIdx-1,:); 
         x_loc'            zeros(1,nc-1); 
         zeros(nc-selfIdx,n) myEye(selfIdx:end,:)];                                 
                                 
end