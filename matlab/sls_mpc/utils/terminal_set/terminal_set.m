function [params, stats] = terminal_set(sys, params)
% Distributed computation of terminal set (following Alg 10.1 from Borelli)
% sys    : LTISystem (we will use sys.A and sys.B2)
% params : MPCParams object containing locality, constraints, etc.
% 
% Currently this function makes the following assumptions:
% 1. Each state is constrained by an upper and lower bound
% 2. No inputs are constrained
% 3. No coupling induced by state constraints

Nx = sys.Nx;

commsAdj  = abs(sys.A) > 0;
stateSupp = commsAdj^(params.locality_-1) > 0;
s_c   = get_locality_col(stateSupp);
s_r   = get_locality_row(stateSupp);

% Initialize with state constraints
HTerm = [params.stateConsMtx_; -params.stateConsMtx_]; 
hTerm = [params.stateUB_; -params.stateLB_];

% Runtime and iterations
times = zeros(Nx, 1); % per state
iters = 0;

% Unconstrained infinite-horizon Phi for precursor set calculation
PhiH2 = h2_sls(sys, params);

Eye  = eye(Nx);
Phi1 = zeros(Nx);

% TODO: there is probably a way to simplify this
for i = 1:Nx % Distributed across columns
    tic;
    Phi1(s_c{i}, i) = PhiH2{i} * Eye(s_c{i},i);
    times(i) = times(i) + toc;
end  

while true
    iters = iters+1;
    
    % Precursor set of current terminal set defined by HNew*x <= hNew
    HNew = zeros(size(HTerm));
    hNew = zeros(size(hTerm));

    for i = 1:Nx
        tic;
        % Two columns: corresponding to upper and lower bounds
        
        % TODO: possibly don't have to do full multiplication
        % i.e. can use Phi1(s_c{i}, i)??
        HNew(s_c{i}, i)    = HTerm(s_c{i}, :)*Phi1(:,i);
        HNew(s_c{i}+Nx, i) = HTerm(s_c{i}+Nx, :)*Phi1(:,i);
        
        hNew(i)    = hTerm(i);
        hNew(i+Nx) = hTerm(i+Nx);
        times(i)   = times(i) + toc;
    end
     
    % Check which rows of HNew are redundant compared to HTerm
    rHNew  = assign_rows_h_terminal_only(sys, params, HNew);
    rHTerm = assign_rows_h_terminal_only(sys, params, HTerm);
        
    redundantRows = [];
    for i=1:sys.Nx
        tic;
        HLocTerm = HTerm(rHTerm{i}, s_r{i});
        HLocNew  = HNew(rHNew{i}, s_r{i});
        hLocTerm = hTerm(rHTerm{i});
        hLocNew  = hNew(rHNew{i});
        
        nRowsNew = length(rHNew{i}); % Rows to check for redundancy
        
        % Concatenate the two
        HLoc = [HLocNew; HLocTerm];
        hLoc = [hLocNew; hLocTerm];
        
        redundantRowsLoc = []; % Indexed according to HLoc

        % Always put the row-to-check first; and only check against rows
        % that are yet to be checked or determined to be not-redundant
        for j=1:nRowsNew
            rows = setdiff(1:size(HLoc, 1), [redundantRowsLoc j]);

            if is_redundant(HLoc([j rows], :), hLoc([j rows]), 1)
                redundantRows(end+1)    = rHNew{i}(j);
                redundantRowsLoc(end+1) = j;
            end
        end
        times(i) = times(i) + toc;
    end
        
    % Remove redundant rows + do intersection
    HNew(redundantRows, :) = [];
    hNew(redundantRows)    = [];
 
    HTerm = [HTerm; HNew];
    hTerm = [hTerm; hNew];
    
    if isempty(HNew)
        % Set has converged: all new rows are redundant
        break; % Stop iterating
    end
end

% Running stats (runtime, iters)
stats         = MPCStats();
stats.time_   = mean(times); % average across all states
stats.iters_  = iters;

params.terminal_H_ = HTerm;
params.terminal_h_ = hTerm;

end