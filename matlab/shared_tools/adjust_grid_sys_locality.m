function AComm = adjust_grid_sys_locality(A)
% SLS toolbox functions determine communication structures in the system
%     by assuming each node corresponds to a single state, and preserving
%     directionality of state/input influences from system matrices A, B.
%     This can lead to undercommunication:
%       1. If state/input influences are uni-directional/asymmetric
%       2. If each node has more than one state
% This function provides an adjusted communication matrix A to take the
% place of (A ~= 0). This matrix reflects the communication structure
%     corresponding to assumptions for the power grid model (2 nodes per
%     state, symmetric communication)

AComm    = A ~= 0;
numNodes = length(A) / 2;

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
        if A(freqi, phasej) ~= 0
            % All possible communications between i and j should exist
            AComm(phasei, phasej) = 1;
            AComm(phasei, freqj)  = 1;
            
            % adjA(freqi, phasej) = 1; % This one already exists
            AComm(freqi, freqj)   = 1;
            
            AComm(phasej, phasei) = 1;
            AComm(phasej, freqi)  = 1;
            
            AComm(freqj, phasei)  = 1;
            AComm(freqj, freqi)   = 1;
        end
    end
end
