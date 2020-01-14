function gap = get_constr_gap(sys, cParams, slsOuts, ctrller)
% Gives roughly ||[Rc; Mc] - [R; M][zI-A -B][Rc; Mc]||_frob
% (Slightly adjusted due to Rc{1} being forced to be identity)
% this is the entity that should be minimized by using a 
% least-squares approach

F  = get_ctrller_constraint(sys, slsOuts, cParams.T_);
F1 = F(:, 1:sys.Nx);
F2 = F(:,sys.Nx+1:end);

RcBlock = zeros(sys.Nx*cParams.T_, sys.Nx);
McBlock = zeros(sys.Nu*cParams.T_, sys.Nx);

% below is copied from get_ctrller_constraint.m
for i=1:cParams.T_
    ix = get_range(i, sys.Nx);
    iu = get_range(i, sys.Nu);        
    RcBlock(ix, :) = ctrller.Rc_{i};
    McBlock(iu, :) = ctrller.Mc_{i};
end

RMcBlock = [RcBlock(sys.Nx+1:end, :); McBlock];

gap = norm(F2 * RMcBlock + F1, 'fro');
end
