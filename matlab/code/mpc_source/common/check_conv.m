criterion_failed = false;

for sys = 1:Nx
      local_phi = Phi(r{sys},s_r{sys}(tFIR,:));
      local_psi = Psi(r{sys},s_r{sys}(tFIR,:));
      local_psi_prev = Psi_prev(r{sys},s_r{sys}(tFIR,:));

      local_diff_d = norm(local_psi-local_psi_prev,'fro');
      local_diff_p = norm(local_phi-local_psi,'fro');
            
      if local_diff_p > eps_p || local_diff_d > eps_d
          criterion_failed = true;
          break; % if one fails, can stop checking the rest
      end
end
