function M = make_m(params, Phi_x, Phi_u)

Nx = size(Phi_x{1}, 1); T = length(Phi_x); % Doesn't use tFIR from params

M = zeros(Nx, Nx);
for k=1:T
    M = M + abs(params.Hx_*Phi_x{k} + params.Hu_*Phi_u{k}); % element-wise abs
end

end