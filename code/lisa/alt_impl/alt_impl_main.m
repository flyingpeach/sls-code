%% choose the system you want to work with
setup1;

% for rank/zero conditions, try to match the precision of cvx_precision low
% http://cvxr.com/cvx/doc/solver.html#solver-precision
eps = 2.22e-16;
tol = eps.^(3/8);

close all;
%% sandbox
Tc = 15;
%slsOuts_alt = find_alt_impl_block(sys, slsParams, slsOuts, Tc);
slsOuts_alt = find_alt_impl_precise(sys, slsParams, slsOuts, Tc);

s_a{1}  = slsOuts_alt;

met     = AltImplMetrics(Tc, tol);
met     = calc_mtx_metrics(met, sys, slsParams, slsOuts, s_a);
met     = calc_cl_metrics(met, sys, simParams, slsParams, slsOuts, s_a);
met

visualize_matrices(slsOuts, slsOuts_alt, Tc, 'all');

% check solver/feasibility statuses
disp(['Statuses:', print_statuses(sys, slsParams, slsOuts, s_a, tol)]);

%% find new implementations Rc, Mc
Tcs    = 2:slsParams.tFIR_+5;
numTcs = length(Tcs);

% slsOuts_alts{i} contains alternate implementation for Tc=Tcs(i)
slsOuts_alts = cell(numTcs, 1);

for idx=1:length(Tcs)
    Tc = Tcs(idx);
    %slsOuts_alts{idx} = find_alt_impl_block(sys, slsParams, slsOuts, Tc);
    slsOuts_alts{idx} = find_alt_impl_precise(sys, slsParams, slsOuts, Tc);
end

% check feasibility / solver statuses
disp(['Statuses:', print_statuses(sys, slsParams, slsOuts, slsOuts_alts, tol)]); 

%% plot stuff
% we might not want to plot all of slsOuts_alts, so this is the sandbox to
% adjust slsOuts_alts / Tcs to our liking

TcsWanted = Tcs;
myIdx = [];

for i=1:numTcs
    Tc = Tcs(i);
    if ismember(Tc, TcsWanted)
        myIdx = [myIdx, i];
    end
end

Tcs_plot = Tcs(myIdx);
slsOuts_alts_plot = slsOuts_alts(myIdx);

met = AltImplMetrics(Tcs_plot, tol);

% calculate matrix-specific metrics
met = calc_mtx_metrics(met, sys, slsParams, slsOuts, slsOuts_alts_plot);

% calculate closed-loop metrics and plot select heat maps
met = calc_cl_metrics(met, sys, simParams, slsParams, slsOuts, slsOuts_alts_plot);

% plot metrics and save to file
savepath = 'C:\Users\Lisa\Desktop\caltech\research\implspace\tmp\';
plot_metrics(met, savepath);

% print again
disp(['Statuses:', print_statuses(sys, slsParams, slsOuts, slsOuts_alts, tol)]);  