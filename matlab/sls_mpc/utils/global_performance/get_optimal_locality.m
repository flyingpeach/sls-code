function [locality, parTime, rankTime, rankDef] = get_optimal_locality(sys, params, varargin)
% locality: size of local communication region (d-hop neighbors) for which
%           the size of the trajectory space is unchanged compared to 
%           global communication, assuming x0 is dense
%           convention: d=1 means self-communication only
% parTime : time (seconds) taken for parallel part of algorithm
% rankTime: time (seconds) taken for rank determination (not parallel right now)
% rankDef : whether we encountered rank deficiency

% params : MPCParams(); the locality_ field will be populated with the 
%          ideal locality
% optional input eps: used to determine whether solution exists for local
%                     subspace; default value is 1e-8


if ~isempty(varargin) > 0
    eps = varargin{1}; 
else
    eps = 1e-8; % Default value
end

% Values of x0 need to be nonzero; specific value doesn't matter
x0 = ones(sys.Nx, 1);

parTime  = 0;
rankTime = 0;
rankDef  = false;

maxLoc   = sys.Nx;
for locality=2:maxLoc
    fprintf('Checking locality size %d\n', locality);
    params.locality_ = locality;
    
    [mtx, parTime1] = get_local_subspace(sys, x0, params, eps);
    tic;
    rankRatio = rank(full(mtx)) / (sys.Nu*(params.tFIR_-1));
    rankTime1 = toc;
    
    parTime  = parTime + parTime1;
    rankTime = rankTime + rankTime1;
    
    if rankRatio == 1
        break;
    elseif rankRatio > 0
        fprintf('Rank deficiency encountered, rank ratio = %d\n', rankRatio); % This is a rare case, give printout
        rankDef = true;
    end
end