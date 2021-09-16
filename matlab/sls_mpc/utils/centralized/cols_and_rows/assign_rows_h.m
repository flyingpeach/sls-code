function r = assign_rows_h(sys, params)
% r{i} represents the set of rows subsystem i solves for

Nx = sys.Nx; Nu = sys.Nu; T = params.tFIR_;
statesUB = []; % which states have an upper bound defined
inputsUB = [];
statesLB = [];
inputsLB = [];

for i = 1:Nx
    if params.has_state_cons() && ~all(params.stateConsMtx_(i,:) == 0)
        if ~isinf(params.stateUB_(i))
            statesUB = [statesUB i];
        end
        if ~isinf(params.stateLB_(i))
            statesLB = [statesLB i];
        end
    end
end

for i = 1:Nu
    if params.has_input_cons() && ~all(params.inputConsMtx_(i,:) == 0)
        if ~isinf(params.inputUB_(i))
            inputsUB = [inputsUB i];
        end
        if ~isinf(params.inputLB_(i))
            inputsLB = [inputsLB i];
        end
    end
end

r = cell(Nx, 1);
nStatesUB = length(statesUB); % H is stacked in this order
nInputsUB = length(inputsUB);
nStatesLB = length(statesLB);
nInputsLB = length(inputsLB);

for i = 1:Nx
    actuator = find(sys.B2(i,:)); % assumes 1 actuator per system    

    if ismember(i, statesUB) 
        for t=1:T
            r{i}(end+1) = nStatesUB*(t-1) + find(statesUB == i);
        end
    end
    
    rowOffset = nStatesUB*T;
    
    if ~isempty(actuator) % add rows for actuator too
        if ismember(actuator, inputsUB)
            for t=1:T-1
                r{i}(end+1) = rowOffset + nInputsUB*(t-1) + find(inputsUB == actuator);
            end
        end
    end
    
    rowOffset = nStatesUB*T + nInputsUB*(T-1);
            
    if ismember(i, statesLB)
        for t=1:T
            r{i}(end+1) = rowOffset + nStatesLB*(t-1) + find(statesLB == i);
        end
    end
    
    rowOffset = nStatesUB*T + nInputsUB*(T-1) + nStatesLB*T;

    if ~isempty(actuator) % add rows for actuator too
        if ismember(actuator, inputsLB)
            for t=1:T-1
                r{i}(end+1) = rowOffset + nInputsLB*(t-1) + find(inputsLB == actuator);
            end
        end
    end
end

end