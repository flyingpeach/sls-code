function intStabMtx = get_int_stab_mtx(sys, ctrller)

Rc = ctrller.Rc_; 
Mc = ctrller.Mc_;
Tc = length(Rc);

Dellc = cell(Tc, 1);

for k=1:Tc-1
    Dellc{k} = Rc{k+1} - sys.A*Rc{k} - sys.B2*Mc{k};
end
Dellc{Tc} = - sys.A*Rc{Tc} - sys.B2*Mc{Tc};

intStabMtx = zeros(sys.Nx*Tc, sys.Nx*Tc);

for i=1:Tc-1
    ix     = get_range(i, sys.Nx);
    ixplus = get_range(i+1, sys.Nx); 
    intStabMtx(ix, ixplus) = eye(sys.Nx);
end

endx = get_range(Tc, sys.Nx);
for i=1:Tc
    ix   = get_range(i, sys.Nx);    
    intStabMtx(endx, ix) = -Dellc{Tc+1-i};
end
