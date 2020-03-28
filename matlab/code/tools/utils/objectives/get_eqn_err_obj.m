function objective = get_eqn_err_obj(sys, clMaps, Rc, Mc)
% Objective for equation error (heuristic difference between desired 
% closed-loop maps and implemented closed-loop maps)
Tc = length(Rc);

[F, G] = get_ctrller_eqn(sys, clMaps, Tc);

RcMc = [];

for k=2:Tc % Omit Rc{1} (assumed to be identity)
    RcMc = [RcMc; Rc{k}];
end
for k=1:Tc
    RcMc = [RcMc; Mc{k}];
end

eqnErr = F*RcMc - G;
objective = norm(eqnErr, 'fro');

end