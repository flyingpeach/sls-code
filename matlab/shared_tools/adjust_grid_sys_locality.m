function adjSys = adjust_grid_sys_locality(sys)
% SLS toolbox functions determine communication structures in the system
%     by assuming each node corresponds to a single state, and preserving
%     directionality of state/input influences from system matrices A, B.
%     This can lead to undercommunication:
%       1. If state/input influences are uni-directional/asymmetric
%       2. If each node has more than one state
% This function provides an adjusted system with an adjusted matrix A,
%     which, when treated by automatic communication structure detection
%     functions, will give out the appropriate communication structure
%     corresponding to assumptions for the power grid model (2 nodes per
%     state, symmetric communication)

adjSys   = sys.copy(); % deepcopy
numNodes = sys.Nx/2;

CONST = -1; % What to populate extra entries with

for i=1:numNodes
    for j=1:numNodes
        if i==j
            continue;
        end
        phasei = 2*i - 1;
        freqi  = 2*i;
        phasej = 2*j - 1;
        freqj  = 2*j;
        
        % This is specific to swing-equation grid system
        if sys.A(freqi, phasej) ~= 0
            % All possible communications between i and j should exist
            adjSys.A(phasei, phasej) = CONST;
            adjSys.A(phasei, freqj)  = CONST;
            
            % adjSys.A(freqi, phasej) = CONST; % This one already exists
            adjSys.A(freqi, freqj)   = CONST;
            
            adjSys.A(phasej, phasei) = CONST;
            adjSys.A(phasej, freqi)  = CONST;
            
            adjSys.A(freqj, phasei)  = CONST;
            adjSys.A(freqj, freqi)   = CONST;
        end
    end
end
