function w = get_loc_bnd_noise(w, sys, params, sigma)
% Truncates given disturbance w into locally bounded disturbance
% Locally bounded disturbance is bounded || ||_2 <= sigma

d        = params.locality_;
T        = params.tFIR_;
tHorizon = size(w, 2);

commsAdj   = abs(sys.A) > 0;
RSupp      = commsAdj^(d-1) > 0;

% Maximum number of neighbors
maxLoc = 0;
for i=1:sys.Nx
    locIdx = union(find(RSupp(i,:)), find(RSupp(:,i)));
    maxLoc = max(maxLoc, length(locIdx));
end

% Maximum norm per element to achieve total patch norm of sigma
% Since frobenius norm squared is additive 
maxNorm = sigma / sqrt(maxLoc) / sqrt(T-1);
for t=1:tHorizon
    for i=1:sys.Nx
        if abs(w(i,t)) > maxNorm
            w(i,t) = w(i,t) / abs(w(i,t)) * maxNorm;
        end
    end
end

% Sanity check
for i=1:sys.Nx
    locIdx = union(find(RSupp(i,:)), find(RSupp(:,i)));
    for t=1:tHorizon-T+2
        tIdx = t:t+T-2; % width = T-1
        if norm(w(locIdx, tIdx), 'fro') > sigma
            error('Noise is not locally bounded!');
        end
    end
end

end