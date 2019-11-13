function plot_ctrller_matrices(sys, slsOuts, ctrller, t)
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

R  = slsOuts.R_;  M  = slsOuts.M_;
Rc = ctrller.Rc_; Mc = ctrller.Mc_;

if strcmp(t, 'all')
    T  = length(slsOuts.R_);
    Tc = length(ctrller.Rc_);

    ts = 1:max(T, Tc);
    statusTxt = sprintf('T=%d, Tc=%d', T, Tc);
    disp(statusTxt);
    
    % zero pad for convenient plotting comparison
    for t_=T+1:Tc
        R{t_} = zeros(sys.Nx, sys.Nx);
        M{t_} = zeros(sys.Nu, sys.Nx);
    end
    for t_=Tc+1:T
        Rc{t_} = zeros(sys.Nx, sys.Nx);
        Mc{t_} = zeros(sys.Nu, sys.Nx);
    end

else
    ts = t;
end

figure; 

for t_=ts % either a single value or array
    subplot('position', pos11);
    imagesc(log10(abs(Rc{t_}))); 
    colorbar; colormap jet; axis equal; axis tight;
    title(sprintf('log10(|Rc|) k=%d', t_)); caxis([logmin logmax]);

    subplot('position', pos12);
    imagesc(log10(abs(R{t_}))); 
    colorbar; colormap jet; axis equal; axis tight;
    title(sprintf('log10(|R|) k=%d', t_)); caxis([logmin logmax]);

    subplot('position', pos21);
    imagesc(log10(abs(sys.B2*Mc{t_}))); 
    colorbar; colormap jet; axis equal; axis tight;
    title(sprintf('log10(|B*Mc|) k=%d', t_)); caxis([logmin logmax]);

    subplot('position', pos22);
    imagesc(log10(abs(sys.B2*M{t_}))); 
    colorbar; colormap jet; axis equal; axis tight;
    title(sprintf('log10(|B*M|) k=%d', t_)); caxis([logmin logmax]);

    pause(0.7);
end
end
