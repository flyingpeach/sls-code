function objective = traj_track(R, M, TFIR, Nx, Nu, xDesired)
% reformulation of (x-xd)' Q (x-xd) + u' R u
% i.e. track desired trajectory while minimizing input

% TODO: should be able to pass these in
xPenalty = eye(Nx); % matrix Q in cost function
uPenalty = 0.1 * eye(Nu); % matrix R in cost function

% TODO: temp hack for proof-of-concept; will need to max over w
w = ones(Nx, 40);

objective = 0;
for t = 1:TFIR
    % need to do the vect operation because of quirk in cvx
    vect = vec(blkdiag(sqrtm(xPenalty), sqrtm(uPenalty)) * ([R{t}; M{t}] * w(:,t) - [xDesired(:,t); zeros(Nu, 1)]));

    objective = objective + vect'*vect;
end
