function sanity_check_terminal_coupling(sys, d, H)
% Terminal set is defined by H*x <= h
% Coupling is only allowed to occur within local sets (which is specified
% by the locality constraint); check that this is the case
% Note: This only checks state coupling, not input coupling

Comms_Adj  = abs(sys.A) > 0;
% Sparsity patterns given by locality constraint d
StateSupp  = Comms_Adj^(d-1) > 0;

for i=1:size(H, 1)
    node = find(H(i, :), 1); % first nonzero entry
    
    % Note: this assumes symmetry in A
    % Nonzero indices of H_term should be subset of StateSupp
    Comms_Hi = H(i, :) ~= 0; % boolean
    if ~isempty(find(StateSupp(node, :) - Comms_Hi < 0, 1))
        mpc_error('Terminal set coupling violates locality!');
    end
end

end