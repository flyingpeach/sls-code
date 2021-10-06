function [H_term, h_term, stats] = terminal_set(sys, params)
% Distributed computation of terminal set
% sys    : LTISystem (we will use sys.A and sys.B2)
% params : MPCParams object containing locality, constraints, etc.
% 
% Currently this function makes the following assumptions:
% 1. Each state is constrained by an upper and lower bound
% 2. No inputs are constrained
% 3. No coupling induced by state constraints

Nx = sys.Nx; A = sys.A;
d = params.locality_;

sparsity_x    = abs(A^(d-1)) > 0;    
row_neighbors = cell(Nx, 1);
col_neighbors = cell(Nx, 1);
for i = 1:Nx
    row_neighbors{i} = find(sparsity_x(i, :));
    col_neighbors{i} = find(sparsity_x(:, i));
end

Phi_x  = eye(Nx);
H_term = []; 
h_term = [];
redundant_rows = [];

H = [params.stateConsMtx_; -params.stateConsMtx_];
h = [params.stateUB_; -params.stateLB_];

% Compute unconstrained infinite-horizon Phi
Phi_x_H2 = h2_sls(sys, params);

% Runtime and iterations
times = zeros(Nx, 1); % per state
iters = 0;

% Iterate until the set converges, i.e. all new rows are redundant
while length(redundant_rows) < size(H, 1) 
    iters = iters+1;
    
    % k-step precursor set {x: H_new*x <= h_new} to {x: H*x <= h} under optimal controller
    % At the first step, H_new = H (because Phi_x = I)
    H_new = zeros(size(H));
    h_new = zeros(size(h));
    
    for i = 1:Nx
        tic;
        % Two rows corresponding to upper and lower bounds
        H_new(col_neighbors{i},i)    = H(col_neighbors{i},:)*Phi_x(:,i);
        H_new(col_neighbors{i}+Nx,i) = H(col_neighbors{i}+Nx,:)*Phi_x(:,i);
        
        h_new(i)    = h(i);
        h_new(i+Nx) = h(i+Nx);
        times(i)    = times(i) + toc;
    end
    
    % Check which rows of H_new are redundant compared to H

    % Note: ideally we would check H_new compared to H_term
    % What we do here may take longer to converge and give more redundancy
    redundant_rows = [];
    for i = 1:Nx
        tic;
        for j = [i, i+Nx] % Rows corresponding to upper and lower bounds
            
            % Find worst possible x for this local patch based on H
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
                % If x satisfies H*x <= h, it automatically satisfies
                % H_new(j,:)*x <= h_new(j); so this row in H_new is redundant                
                redundant_rows(end+1) = j;
            end
        end
        times(i) = times(i) + toc;
    end
    
    % Remove desired rows
    H_new(redundant_rows, :) = []; 
    h_new(redundant_rows)    = [];
    
    % Set intersection (minus redundant rows)
    H_term = [H_term; H_new];
    h_term = [h_term; h_new];

    % Compute the infinite-horizon closed-loop map recursively
    % according to eq (13) in https://arxiv.org/pdf/2010.02440
    for i = 1:Nx % Distributed across columns
        tic;
        Phi_x(col_neighbors{i},i) = Phi_x_H2{i} * Phi_x(col_neighbors{i},i);
        times(i) = times(i) + toc;
    end  
end

% Running stats (runtime, iters)
stats         = MPCStats();
stats.time_   = mean(times); % average across all states
stats.iters_  = iters;

end