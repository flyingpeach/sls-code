%% choose the system you want to work with
setup1;

%% sandbox
Tc = 20;
slsOuts_alt = find_alt_impl(sys, slsParams, slsOuts, Tc);

s_a{1}  = slsOuts_alt;
zThresh = 1e-6;
met     = AltImplMetrics([Tc], zThresh);
met     = calc_mtx_metrics(met, sys, slsParams, slsOuts, s_a);
met     = calc_cl_metrics(met, sys, simParams, slsParams, slsOuts, s_a);
met

visualize_matrices(slsOuts, slsOuts_alt, Tc, 'all');

%% find new implementations Rc, Mc
Tcs = [2:20];

% slsOuts_alts{i} contains alternate implementation for Tc=Tcs(i)
slsOuts_alts = cell(length(Tcs), 1);

for i=1:length(Tcs)
    Tc              = Tcs(i);
    slsOuts_alts{i} = find_alt_impl(sys, slsParams, slsOuts, Tc, 'approx');
end

disp(['Statuses:', print_statuses(slsOuts_alts)]); % check solver statuses

%% calculate stats & generate plots
zThresh = 1e-6;
met     = AltImplMetrics(Tcs, zThresh);

% calculate matrix-specific metrics
met = calc_mtx_metrics(met, sys, slsParams, slsOuts, slsOuts_alts);

% calculate closed-loop metrics and plot select heat maps
met = calc_cl_metrics(met, sys, simParams, slsParams, slsOuts, slsOuts_alts);

% plot metrics and save to file
savepath = 'C:\Users\Lisa\Desktop\caltech\research\implspace\tmp\';
plot_metrics(met, savepath);