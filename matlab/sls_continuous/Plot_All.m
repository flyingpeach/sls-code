%clc,clear,close all

fs = 12;
figure(1);

% Before running this file, run the appropriate simulations files (e.g.
% Simple_Poles_H2.m for simple pole values associated with H2 control) to
% obtain the relevant data. You can directly run this file afterwards, or
% save the relevant data into a .mat file and load it.

% H2 part
%load("Grid_Model_Costs.mat");
%load('H2w_repeatedPoles.mat');
%load('H2w_repeatedPoles_2.mat');
SPA_Ratio = SPA_SLS_Ratio;
RP_Ratio = Repeated_Poles_Ratios;
RP_Ratio2 = Repeated_Poles_Ratios2;
Sparse_Ratio = Sparse_SLS_Ratio;

subplot(2,3,1); axis square; box on;
hold on 
plot(Mult,RP_Ratio2,'gx-',LineWidth=1.25);
plot(PoleNumber(1:4),SPA_Ratio(1:4),"bo-",LineWidth=2);
plot(Mult,RP_Ratio,'r*-',LineWidth=1.25);
plot(PoleNumber(1:4),ones(1,4),"k:", LineWidth=2);
legend('Ritz, $p=-0.3$', 'SPA', 'Ritz, $p=-1$', 'interpreter','latex')
title('Dense, varying poles','FontSize',12)
ylabel('normalized cost','FontSize',12)
ax = gca;
ax.XAxis.FontSize = fs;
ax.YAxis.FontSize = fs;
xlim([1.6 8.4])
title({'cost vs. # of poles', 'centralized'},'FontSize',12);
ylim([0.7 5.2])

subplot(2,3,2); axis square; box on;
hold on
plot(PoleNumber(1:4),Sparse_SPA_Ratio(1:4),"bo-",LineWidth=2);
plot(PoleNumber(1:4),ones(1,4),"k:", LineWidth=2);
title('d=2, varying poles','FontSize',12)
ax = gca;
ax.XAxis.FontSize = fs;
ax.YAxis.FontSize = fs;
xlim([1.6 8.4])
title({'cost vs. # of poles', 'd=2'},'FontSize',12);
ylim([0.7 5.2])

subplot(2,3,3); axis square; box on;
hold on
plot(Dist(1:5),Sparse_Ratio(1:5),"bo-",LineWidth=2);
plot(Dist(1:5),ones(1,5),"k:", LineWidth=2);
title('Fixed poles, varying d','FontSize',12)
ax = gca;
ax.XAxis.FontSize = fs;
ax.YAxis.FontSize = fs;
title({'cost vs. d', '# of poles = 6'},'FontSize',12);
xlim([0.78 5.22])
ylim([0.7 5.2])

% HInf part
load("Grid_Model_Hinf_Costs.mat");
SPA_Ratio = SPA_SLS_Ratio;
Sparse_Ratio = Sparse_SLS_Ratio;

subplot(2,3,4); axis square; box on;
hold on 
plot(PoleNumber,SPA_Ratio,"bo-",LineWidth=2.0);
plot(PoleNumber,ones(1,size(PoleNumber,1)),"k:", LineWidth=2.0);

xlabel('# of poles','FontSize',12)
ylabel('normalized cost','FontSize',12)
ax = gca;
ax.XAxis.FontSize = fs;
ax.YAxis.FontSize = fs;
xlim([1.6 8.4])
ylim([0.95 1.75])

subplot(2,3,5); axis square; box on;
hold on
plot(PoleNumber,Sparse_SPA_Ratio,"bo-",LineWidth=2.0);
plot(PoleNumber,ones(1,size(PoleNumber,1)),"k:", LineWidth=2.0);
xlabel('# of poles','FontSize',12)
ax = gca;
ax.XAxis.FontSize = fs;
ax.YAxis.FontSize = fs;
xlim([1.6 8.4])
ylim([0.95 1.75])

subplot(2,3,6); axis square; box on;
hold on
plot(Dist(1:5),Sparse_Ratio(1:5),"bo-",LineWidth=2.0);
plot(Dist(1:5),ones(1,5),"k:", LineWidth=2.0);
title({'cost vs. d', '# of poles = 4'},'FontSize',12);
xlabel('d','FontSize',12)
ax = gca;
ax.XAxis.FontSize = fs;
ax.YAxis.FontSize = fs;
xlim([0.78 5.22])
ylim([0.95 1.75])