function sanity_check_actuation(sys)
% Implementation of distributed MPC currently requires
%  1. One actuator per subsystem
%  2. Actuator indices should be a strictly increasing
%     function of the subsystem indices they correspond to
%  Example: B2 = [0 0 1;
%                [0 1 0] violates (2)
%           since actuator 3 corresponds to subsystem 1
%           but   actuator 2 corresponds to subsystem 2
%  B2 = [0 1 0;
%       [0 0 1] will work


lastActIdx = 0;
for row=1:sys.Nu
    if length(find(sys.B2(row,:))) > 1
        mpc_error('Maximum one actuator per subsystem');
    end
    if find(sys.B2(row,:)) <= lastActIdx
        mpc_error('sys.B2 has unsupported sparsity pattern');
    end
    lastActIdx = find(sys.B2(row,:));
end

end