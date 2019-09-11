function visualize_matrices(slsOuts, slsOuts_alt, Tc, t)
% Plot log magnitudes of R, M, Rc, Mc
% Inputs
%     slsOuts     : original closed loop system (contains R, M)
%     slsOuts_alt : alternate CL implementations (contains Rc, Mc)
%     Tc          : the Tc for which slsOuts_alt was calculated
%     t           : which timestep to plot R/M/Rc/Mc for
%                   set 'all' if you want an animated scroll-through

logmin = -6; logmax = 1;

Nu      = size(slsOuts.M_{1}, 1);
Nx      = size(slsOuts.M_{1}, 2);

% nice spacing / arrangement
xgap    = 0.05;
yBuffer = 0.02;

xfill  = 1 - 4 * xgap;
sizeNx = xfill / 2;
sizeNu = sizeNx * Nu / Nx;
yspace = max(1 - sizeNx - sizeNu - yBuffer, 0);
ygap   = yspace / 4;

xpos1 = xgap; xpos2 = sizeNx + 3 * xgap;
ypos1 = ygap; ypos2 = sizeNu + 3 * ygap;

% [left bottom width height]
pos11 = [xpos1 ypos1 sizeNx sizeNu];
pos21 = [xpos2 ypos1 sizeNx sizeNu];
pos12 = [xpos1 ypos2 sizeNx sizeNx];
pos22 = [xpos2 ypos2 sizeNx sizeNx];

if strcmp(t, 'all')
    T = length(slsOuts.R_);
    ts = 1:min(T, Tc);
    statusTxt = sprintf('T=%d, Tc=%d; animating for %d timesteps', T, Tc, min(T, Tc));
    disp(statusTxt);
else
    ts = t;
end

figure; 

for t_=ts % either a single value or array 
    st = suptitle(sprintf('log10(|val|) @ t=%d, Tc=%d', t_, Tc));
    set(st, 'FontSize', 12);
    
    subplot('position', pos12);
    imagesc(log10(abs(slsOuts_alt.R_{t_}))); 
    colorbar; colormap jet; axis equal; axis tight;
    title('Rc'); caxis([logmin logmax]);

    subplot('position', pos22);
    imagesc(log10(abs(slsOuts.R_{t_}))); 
    colorbar; colormap jet; axis equal; axis tight;
    title('R'); caxis([logmin logmax]);

    subplot('position', pos11);
    imagesc(log10(abs(slsOuts_alt.M_{t_}))); 
    colorbar; colormap jet; axis equal; axis tight;
    title('Mc'); caxis([logmin logmax]);

    subplot('position', pos21);
    imagesc(log10(abs(slsOuts.M_{t_}))); 
    colorbar; colormap jet; axis equal; axis tight;
    title('M'); caxis([logmin logmax]);

    pause(0.7);
end
end
