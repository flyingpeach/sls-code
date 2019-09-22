function spectralRadius = check_int_stability(sys, Tc, Rc, Mc)
for k=1:Tc-1
    Dellc{k} = Rc{k+1} - sys.A*Rc{k} - sys.B2*Mc{k};
end
Dellc{Tc} = - sys.A*Rc{Tc} - sys.B2*Mc{Tc};

intStabMtx = zeros(sys.Nx*(Tc+1), sys.Nx*(Tc+1));

% copied from get_F
get_range = @(idx, size) (size*(idx-1)+1:size*(idx-1)+size);

for i=1:Tc
    ix     = get_range(i, sys.Nx);
    ixplus = get_range(i+1, sys.Nx); 
    intStabMtx(ix, ixplus) = eye(sys.Nx);
end

for i=2:Tc
    ix   = get_range(i, sys.Nx);
    endx = get_range(Tc+1, sys.Nx);
    
    intStabMtx(endx, ix) = Dellc{Tc+2-i};
end

spectralRadius = max(abs(eig(intStabMtx))); % biggest eigenvalue
