function sys = generate_continuous_time_grid_plant(actuatedNodes, adjMtx, susceptMtx, inertiasInv, dampings)
% Given topology of power grid, return system matrices
% Equations are continuous time versions of Anderson et al's power system model (5.10)
%
% actuatedNodes: indices of nodes with actuation
% adjMtx       : adjacency matrix of grid
% susceptMtx   : matrix of susceptance (imaginary admittance) values between nodes
%
% sys          : LTISystem containing system matrices
%

% Only sys.A and sys.B2 are tested for this function. Other system
% parameters may be wrong !!!
numNodes = size(adjMtx, 1);

sys    = LTISystem;
Nx = 2 * numNodes; 
Nu = length(actuatedNodes);
Nw = 2 * numNodes;
Nz = Nx + Nu;

% B1, C1, D11 and D12 are needed when LTISystem is needed
sys.A   = sparse(Nx);
sys.B1  = speye(Nx); % anyone can receive a disturbance
sys.B2  = sparse(Nx, Nu);
sys.C1  = [speye(Nx); sparse(Nu, Nx)];
sys.D11 = sparse(Nz, Nw);
sys.D12 = [sparse(Nx, Nu); speye(Nu)];

actCounter = 1;
for i=1:numNodes
    thisIdx = (2*i-1):(2*i);
    
    mi = inertiasInv(i); % This is actually 1/(m_i)
    di = dampings(i);
    ki = sum(susceptMtx(i,:));
    
    sys.A(thisIdx, thisIdx) = [0 1; -ki*mi -di*mi];
    if ismember(i, actuatedNodes)
        sys.B2(thisIdx, actCounter) = [0; 1]; %ThisIdx is a bunch of indexes
        actCounter = actCounter + 1;
    end
        
    for j=1:numNodes
        if adjMtx(i,j) == 1
            neighbIdx = (2*j-1):(2*j);
            kij = susceptMtx(i,j);
            susceptMtx(i,j) = kij;
            sys.A(thisIdx, neighbIdx) = [0 0; kij*mi 0];
        end
    end
end
end