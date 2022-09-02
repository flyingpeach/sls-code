function PsiSupp = get_sparsity_psi(sys, params)
% This is the sparsity of the entire Toeplitz matrix representing system
% responses. Often only the first block-column (Nx columns) are required.

if ~isempty(sys.AComm)
    Comms_Adj = sys.AComm; % Use the specified communication structure
else
    Comms_Adj = sys.A ~= 0; % Use comm structure implied by A
end

Nx = sys.Nx; Nu = sys.Nu; T = params.tFIR_;   

% Sparsity / support for each block matrix of Psi
RSupp     = Comms_Adj^(params.locality_-1) > 0;
MSupp     = (abs(sys.B2)' * RSupp) > 0;

% Get sparsity of big block matrices
PsixSupp = false(Nx*T, Nx*T);
PsiuSupp = false(Nu*(T-1), Nx*T);
for k = 1:T % block-row of block-toeplitz matrix
    for l = 1:k
        kx = get_range(k, Nx);
        lx = get_range(l, Nx);
        PsixSupp(kx, lx) = RSupp;
        
        if k < T
            ku = get_range(k, Nu);
            lu = get_range(l, Nx);
            PsiuSupp(ku, lu) = MSupp;
        end
    end    
end

PsiSupp = [PsixSupp; PsiuSupp];
end
