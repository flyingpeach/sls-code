function spectralRadius = check_int_stability(sys, Rc, Mc)

Tc = length(Rc);

for k=1:Tc-1
    Dellc{k} = Rc{k+1} - sys.A*Rc{k} - sys.B2*Mc{k};
end
Dellc{Tc} = - sys.A*Rc{Tc} - sys.B2*Mc{Tc};

intStabMtx = zeros(sys.Nx*Tc, sys.Nx*Tc);

% copied from get_F
get_range = @(idx, size) (size*(idx-1)+1:size*(idx-1)+size);

for i=1:Tc-1
    ix     = get_range(i, sys.Nx);
    ixplus = get_range(i+1, sys.Nx); 
    intStabMtx(ix, ixplus) = eye(sys.Nx);
end

for i=1:Tc
    ix   = get_range(i, sys.Nx);
    endx = get_range(Tc, sys.Nx);
    
    intStabMtx(endx, ix) = Dellc{Tc+1-i};
end
spectralRadius = max(abs(eig(intStabMtx))); % biggest eigenvalue
