function [H_term, h_term] = terminal_set(sys, params)
% Distributed computation of terminal set
% sys    : LTISystem (we will use sys.A and sys.B2)
% params : MPCParams object containing locality, constraints, etc.
% 
% Currently this function makes the following assumptions:
% 1. Each state is constrained by an upper and lower bound
% 2. No inputs are constrained
% 3. No coupling induced by state constraints

Nx = sys.Nx; A = sys.A; d = params.locality_;

% Compute unconstrained infinite-horizon Phi
[Phi_x_H2, col_neighbors] = h2_sls(sys, params);

sparsity_x    = abs(A^(d-1)) > 0;
row_neighbors = cell(Nx, 1);
for j = 1:Nx
    row_neighbors{j} = find(sparsity_x(j,:));
end

Phi_x  = eye(Nx);
H_term = []; 
h_term = [];
remove_rows = [];

H = [params.stateConsMtx_; -params.stateConsMtx_];
h = [params.stateUB_; -params.stateLB_];

% Iterate until the set converges, i.e. all new rows are redundant
while length(remove_rows) < size(H, 1) 
    
    % New rows
    H_new = zeros(size(H));
    h_new = zeros(size(h));
    
    for i = 1:Nx
        % Calculate two new rows (corresponding to upper and lower bounds)
        H_new(col_neighbors{i},i)    = H(col_neighbors{i},:)*Phi_x(:,i);
        H_new(col_neighbors{i}+Nx,i) = H(col_neighbors{i}+Nx,:)*Phi_x(:,i);

        h_new(i)    = h(i);
        h_new(i+Nx) = h(i+Nx);
    end
    
    % Share information with neighboring subsystems
    remove_rows = [];
    
    for i = 1:Nx
        for j = [i, i+Nx] % corresponds to upper and lower bounds
            
            % Find worst possible x (i.e. take intersection)
            numNeighbors = length(row_neighbors{i});
            x_worst      = zeros(numNeighbors, 1);

            for l = 1:numNeighbors
                k = row_neighbors{i}(l);                
                if H_new(j,k) > 0
                     x_worst(l) = params.stateUB_(k);
                elseif H_new(j,k) <= 0
                     x_worst(l) = params.stateLB_(k);
                end
            end
            
            if H_new(j, row_neighbors{i})*x_worst < h_new(j)
                % This constraint is redundant; remove it
                remove_rows = [remove_rows, j];
            end
            
        end
    end
    
    % Remove desired rows
    H_new(remove_rows, :) = []; 
    h_new(remove_rows)    = [];
    
    % Intersect the two sets
    H_term = [H_term; H_new];
    h_term = [h_term; h_new];
        
    % Compute the infinite-horizon closed-loop map recursively
    % according to eq (13) in https://arxiv.org/pdf/2010.02440
    for i = 1:Nx % We distribute across columns
        Phi_x(col_neighbors{i},i) = Phi_x_H2{i} * Phi_x(col_neighbors{i},i);
    end
end
