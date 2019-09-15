%% choose the system you want to work with
setup1;

% for rank/zero conditions, try to match the precision of cvx_precision low
% http://cvxr.com/cvx/doc/solver.html#solver-precision
eps = 2.22e-16;
tol = eps.^(1/4);

close all;
%% sandbox
Tc = 2;
slsOuts_alt = find_alt_impl(sys, slsParams, slsOuts, Tc);

s_a{1}  = slsOuts_alt;

met     = AltImplMetrics([Tc], tol);
met     = calc_mtx_metrics(met, sys, slsParams, slsOuts, s_a);
met     = calc_cl_metrics(met, sys, simParams, slsParams, slsOuts, s_a);
met

[rankF, rankF2] = get_rank(sys, slsParams, slsOuts, Tc, tol);
visualize_matrices(slsOuts, slsOuts_alt, Tc, 'all');

% check solver/feasibility statuses
disp(['Statuses:', print_statuses(sys, s_a, rankF, rankF2)]); 

%% find new implementations Rc, Mc
Tcs    = [2:25];
numTcs = length(Tcs);

% slsOuts_alts{i} contains alternate implementation for Tc=Tcs(i)
slsOuts_alts = cell(numTcs, 1);
rankFs       = zeros(numTcs, 1);
rankF2s      = zeros(numTcs, 1);

for idx=1:length(Tcs)
    Tc = Tcs(idx);

    [rankFs(idx), rankF2s(idx)] = get_rank(sys, slsParams, slsOuts, Tc, tol);
    slsOuts_alts{idx} = find_alt_impl(sys, slsParams, slsOuts, Tc);
end

% check solver/feasibility statuses
disp(['Statuses:', print_statuses(sys, slsOuts_alts, rankFs, rankF2s)]); 

%% calculate stats & generate plots
met     = AltImplMetrics(Tcs, tol);

% calculate matrix-specific metrics
met = calc_mtx_metrics(met, sys, slsParams, slsOuts, slsOuts_alts);

% calculate closed-loop metrics and plot select heat maps
met = calc_cl_metrics(met, sys, simParams, slsParams, slsOuts, slsOuts_alts);

% plot metrics and save to file
savepath = 'C:\Users\Lisa\Desktop\caltech\research\implspace\tmp\';
plot_metrics(met, savepath);