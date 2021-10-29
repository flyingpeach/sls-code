function sys = generate_grid_plant(actuatedNodes, adjMtx, susceptMtx, inertiasInv, dampings, Ts)
% Given topology of power grid, return system matrices
% Equations taken from Anderson et al's power system model (5.10)
%
% actuatedNodes: indices of nodes with actuation
% adjMtx       : adjacency matrix of grid
% susceptMtx   : matrix of susceptance (imaginary admittance) values between nodes
% Ts           : sampling time
%
% sys          : LTISystem containing system matrices
%

numNodes = size(adjMtx, 1);

sys    = LTISystem;
sys.Nx = 2 * numNodes; 
sys.Nu = length(actuatedNodes);
sys.Nw = 2 * numNodes;
sys.Nz = sys.Nx + sys.Nu;

sys.A   = zeros(sys.Nx, sys.Nx);
sys.B1  = eye(sys.Nx); % anyone can receive a disturbance
sys.B2  = zeros(sys.Nx, sys.Nu);
sys.C1  = [speye(sys.Nx); sparse(sys.Nu, sys.Nx)];
sys.D11 = sparse(sys.Nz, sys.Nw);
sys.D12 = [sparse(sys.Nx, sys.Nu); speye(sys.Nu)];
sys.sanity_check();

actCounter = 1;
for i=1:numNodes
    thisIdx = (2*i-1):(2*i);
    
    mi = inertiasInv(i);
    di = dampings(i);
    ki = sum(susceptMtx(i,:));
    
    sys.A(thisIdx, thisIdx) = [1 Ts; -ki*mi*Ts 1-di*mi*Ts];
    if ismember(i, actuatedNodes)
        sys.B2(thisIdx, actCounter) = [0; 1]; 
        actCounter = actCounter + 1;
    end
        
    for j=1:numNodes
        if adjMtx(i,j) == 1
            neighbIdx = (2*j-1):(2*j);
            kij = susceptMtx(i,j);
            susceptMtx(i,j) = kij;
            sys.A(thisIdx, neighbIdx) = [0 0; kij*mi*Ts 0];
        end
    end
end