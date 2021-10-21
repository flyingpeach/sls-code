function params = remove_redundancy_terminal(sys, params)
% Returns params with redundant rows in params.terminal_H_ and _h_ removed

rH = assign_rows_h_terminal_only(sys, params);
H  = params.terminal_H_;
h  = params.terminal_h_;

HSizeInit = size(H, 1);

commsAdj      = abs(sys.A) > 0;
stateSupp     = commsAdj^(params.locality_-1) > 0;
redundantRows = [];

for i=1:sys.Nx
    locCols = find(stateSupp(i,:));
    
    HLoc             = H(rH{i}, locCols);
    hLoc             = h(rH{i}, :);
    nRowsLoc         = length(rH{i});
    redundantRowsLoc = []; % indexed according to HLoc

    % Always put the row-to-check first; and only check against rows
    % that are yet to be checked or determined to be not-redundant
    for j=1:nRowsLoc
        rows = setdiff(1:nRowsLoc, [redundantRowsLoc j]);
        
        if is_redundant(HLoc([j rows], :), hLoc([j rows]), 1)
            redundantRows(end+1)    = rH{i}(j);
            redundantRowsLoc(end+1) = j;
        end
    end
end

params.terminal_H_(redundantRows, :) = [];
params.terminal_h_(redundantRows)    = [];

HSizeNow  = size(params.terminal_H_, 1);
HSizeInit = size(H, 1);

fprintf('Removed %d locally redundant constraints\n', HSizeInit-HSizeNow);

end

 