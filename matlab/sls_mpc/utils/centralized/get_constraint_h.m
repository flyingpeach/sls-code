function [H, h] = get_constraint_h(sys, params)
% This function works with time-invariant constraints of the form
%    stateLB <= stateConsMtx*state <= stateUB
%    inputLB <= inputConsMtx*input <= inputUB
% and converts them into the form of H * Phi * [x0; delta] <= h

Nx = sys.Nx; Nu = sys.Nu; T = params.tFIR_;

HStateUB = []; HStateLB = []; hStateUB = []; hStateLB = [];
if params.has_state_cons()
    for t = 1:T
        tx = get_range(t, Nx);
        
        HStateUB(tx, tx) = params.stateConsMtx_;    
        HStateLB(tx, tx) = params.stateConsMtx_;
        hStateUB = [hStateUB; params.stateUB_];
        hStateLB = [hStateLB; params.stateLB_];
    end
    % Make correct number of columns
    HStateUB = [HStateUB zeros(Nx*T, Nu*(T-1))]; 
    HStateLB = [HStateLB zeros(Nx*T, Nu*(T-1))];
end

HInputUB = []; HInputLB = []; hInputUB = []; hInputLB = [];
if params.has_input_cons()
    for t = 1:T-1
        tu = get_range(t, Nu);
        HInputUB(tu, tu+Nx*T) = params.inputConsMtx_;
        HInputLB(tu, tu+Nx*T) = params.inputConsMtx_;
        hInputUB = [hInputUB; params.inputUB_];
        hInputLB = [hInputLB; params.inputLB_];
    end
end

H = [HStateUB; HInputUB; -HStateLB; -HInputLB];
h = [hStateUB; hInputUB; -hStateLB; -hInputLB];

delRow = [];
for i=1:length(h)
    if isinf(h(i)) || all(H(i,:) == 0) % Bound is infinite or constrain coefficient is 0
        delRow = [delRow i]; % Mark that row for removal (trivial constraint)
    end
end
H(delRow, :) = []; % Remove designated rows
h(delRow, :) = [];

end
