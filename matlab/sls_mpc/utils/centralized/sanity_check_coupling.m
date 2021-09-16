function sanity_check_coupling(sys, params)
% Coupling is only allowed to occur within local sets (which is specified
% by the locality constraint); check that this is the case

Comms_Adj  = sys.A ~= 0;
% sparsity patterns given by locality constraint
StateSupp  = Comms_Adj^(params.locality_-1) > 0;
InputSupp  = (abs(sys.B2)' * StateSupp * abs(sys.B2)) > 0;

Comms_Cost_State = params.QSqrt_ ~= 0;
if ~isempty(find(StateSupp - Comms_Cost_State < 0, 1))
    mpc_error('State cost coupling violates locality!');
end

Comms_Cost_Input = params.RSqrt_ ~= 0;
if ~isempty(find(InputSupp - Comms_Cost_Input < 0, 1))
    mpc_error('Input cost coupling violates locality!');
end

if params.has_state_cons()
    Comms_Cons_State = params.stateConsMtx_ ~= 0;
    if ~isempty(find(StateSupp - Comms_Cons_State < 0, 1))
        mpc_error('State constraint coupling violates locality!');
    end 
end

if params.has_input_cons()
    Comms_Cons_Input = params.inputConsMtx_ ~= 0;
    if ~isempty(find(InputSupp - Comms_Cons_Input < 0,1))
        mpc_error('Input constraint coupling violates locality!');
    end
end
        
end