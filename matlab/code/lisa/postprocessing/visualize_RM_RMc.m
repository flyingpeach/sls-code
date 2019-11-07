function visualize_RM_RMc(sys, slsOuts, ctrller, t)
% Plot log magnitudes of R, M, Rc, Mc
% Note: M, Mc are plotted as B2*M and B2*Mc to visualize actuation per node
% Inputs
%     sys     : LTISystem containing system matrices
%     slsOuts : SLSOutputs containing R, M
%     ctrller : Ctrller containing Rc, Mc
%     t       : which timestep to plot R/M/Rc/Mc for
%               set 'all' if you want an animated scroll-through

logmin = -6; logmax = 1;

% nice spacing / arrangement
xgap    = 0.05;
yBuffer = 0.02;

xfill  = 1 - 4 * xgap;
sizeNx = xfill / 2;
yspace = max(1 - 2 * sizeNx - yBuffer, 0);
ygap   = yspace / 4;

xpos1 = xgap; xpos2 = sizeNx + 3 * xgap;
ypos1 = ygap; ypos2 = sizeNx + 3 * ygap;

% [left bottom width height]
pos11 = [xpos1 ypos1 sizeNx sizeNx];
pos21 = [xpos2 ypos1 sizeNx sizeNx];
pos12 = [xpos1 ypos2 sizeNx sizeNx];
pos22 = [xpos2 ypos2 sizeNx sizeNx];

if strcmp(t, 'all')
    T  = length(slsOuts.R_);
    Tc = length(ctrller.Rc_);

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
    imagesc(log10(abs(ctrller.Rc_{t_}))); 
    colorbar; colormap jet; axis equal; axis tight;
    title('Rc'); caxis([logmin logmax]);

    subplot('position', pos22);
    imagesc(log10(abs(slsOuts.R_{t_}))); 
    colorbar; colormap jet; axis equal; axis tight;
    title('R'); caxis([logmin logmax]);

    subplot('position', pos11);
    imagesc(log10(abs(sys.B2*ctrller.Mc_{t_}))); 
    colorbar; colormap jet; axis equal; axis tight;
    title('B*Mc'); caxis([logmin logmax]);

    subplot('position', pos21);
    imagesc(log10(abs(sys.B2*slsOuts.M_{t_}))); 
    colorbar; colormap jet; axis equal; axis tight;
    title('B*M'); caxis([logmin logmax]);

    pause(0.7);
end
end
