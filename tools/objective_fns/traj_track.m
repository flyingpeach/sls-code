function objective = traj_track(sys, params, R, M)
% Reformulation of (x-xd)' Q (x-xd) + u' R u
% i.e. track desired trajectory while minimizing input

% TODO: should be able to pass these in
xPenalty = eye(sys.Nx); % matrix Q in cost function
uPenalty = 0.1 * eye(sys.Nu); % matrix R in cost function

% TODO: temp hack for proof-of-concept; will need to max over w
w = ones(sys.Nx, 40);

objective = 0;
for t = 1:params.tFIR_
    % need to do the vect operation because of quirk in cvx
    vect = vec(blkdiag(sqrtm(xPenalty), sqrtm(uPenalty)) * ([R{t}; M{t}] * w(:,t) - [params.xDes_(:,t); zeros(sys.Nu, 1)]));

    objective = objective + vect'*vect;
end
