close all; clear all; clc;
cd('C:\Users\chri3793\Documents\MATLAB\DPhil\pesgm18_mtlb')
% cd('C:\Users\Matt\Documents\MATLAB\DPhil\pesgm18_mtlb');
addpath('pesgm_funcs');
% fig_loc = [pwd,'\figures\'];
fig_loc = 'C:\Users\chri3793\Documents\MATLAB\DPhil\pesgm18_mtlb\figures\';
% fig_loc = 'C:\Users\Matt\Desktop\wc171023\figures\';

set(0,'defaulttextinterpreter','latex');
set(0,'defaultaxesfontsize',14);
fig_nompossml = [100 100 300 350];
%%
% INPUTS TO DETERMINE BEHAVIOUR -------------------

C = {0.4,[1.0 0.6 0.2],0.4,0.4};
VT = {1.10,1.05,1.00,0.95};

F.pg_ssc = linspace(-1e-6,0.2,200);
F.qg_ssc = linspace(-0.4,1e-6,200);
F.n = 80;
% F.pg_ssc = linspace(-1e-6,0.2,400);
% F.qg_ssc = linspace(-0.4,1e-6,600);
% F.n = 80;

%-------------------------------------------------

% remain the same:
F.SRC = 'SOURCEBUS';
F.filename = '\opendss_models\34Bus\ieee34Mod1_fxd';
F.feeder = '34bus'; % for calc_ww and T_ratios
F.tr_buses = {{'888','832'};{'832','852'};{'850','814'};{'800','SOURCEBUS'}};
F.fig_nompos = [200 150 450 340];

SCD = load('solar_curves_dataset.mat','Ps');
F.Ps0 = SCD.Ps; %pu
F.NUT = '834';
F.Vp = 1.06; %pu

RR = cell(size(C));
for i = 1:numel(C)
    F.V0 = VT{i};
    F.Ps0_k = C{i};
    RR{i} = run_pesgm_feeder( F );
end
%% PLOT #1 + #2: Measured and Estimated Utility versus c
fig = figure('Color','White','Position',fig_nompossml); 
figname = [fig_loc,'msrd_utility'];

plot( RR{2}.Qgd/min(RR{2}.Qgd) , RR{2}.D_Eg ); hold on;
plot( RR{2}.Qgd/min(RR{2}.Qgd) , RR{2}.D_Et ,'--'); hold on;
xlabel('$Q_{g}/\dot{Q}_{g}$'); ylabel('$\Delta E_{(\cdot)}$, UNIT?');
lgnd = legend('$\Delta E_{g}^{\,c\, =\, 1.0}$','$\Delta E_{g}^{\,c\, =\, 0.6}$',...
                '$\Delta E_{g}^{\,c\, =\, 0.2}$','$\Delta E_{t}^{\,c\, =\, 1.0}$',...
                '$\Delta E_{t}^{\,c\, =\, 0.6}$','$\Delta E_{t}^{\,c\, =\, 0.2}$','Location','NorthWest');
set(lgnd,'Interpreter','Latex');
xticks([0 0.25 0.5 0.75 1.0]);
xs = axis; grid on;

% export_fig(fig,figname);
% export_fig(fig,[figname,'.pdf'],'-dpdf');

fig = figure('Color','White','Position',fig_nompossml+[fig_nompossml(3) 0 0 0]); 
figname = [fig_loc,'estd_utility'];

plot( RR{2}.Qgd_est/min(RR{2}.Qgd_est) , RR{2}.D_Eg_est ); hold on;
plot( RR{2}.Qgd_est/min(RR{2}.Qgd_est) , RR{2}.D_Et_est ,'--'); hold on;
xlabel('$Q_{g}/\dot{Q}_{g}$'); ylabel('$\Delta E_{(\cdot)}$, UNIT?');
lgnd = legend('$\Delta E_{g}^{\,c\, =\, 1.0}$','$\Delta E_{g}^{\,c\, =\, 0.6}$',...
                '$\Delta E_{g}^{\,c\, =\, 0.2}$','$\Delta E_{t}^{\,c\, =\, 1.0}$',...
                '$\Delta E_{t}^{\,c\, =\, 0.6}$','$\Delta E_{t}^{\,c\, =\, 0.2}$','Location','NorthWest');
set(lgnd,'Interpreter','Latex'); xticks([0 0.25 0.5 0.75 1.0]);
axis(xs); grid on;

% export_fig(fig,figname);
% export_fig(fig,[figname,'.pdf'],'-dpdf');

%% PLOT #3 + #4: Measured and estimated Error versus c
fig = figure('Color','White','Position',fig_nompossml); 
figname = [fig_loc,'msrd_error'];

plot(RR{2}.Qgd/min(RR{2}.Qgd),RR{2}.e_L); hold on;
xlabel('$Q_{g}/\dot{Q}_{g}$'); ylabel('$\epsilon _{l}, \%$');
xs = axis; grid on;
plot([0.01 0.99],10*[1 1],'k--');
lgnd = legend('$\epsilon _{l}^{\,c\, =\, 1.0}$','$\epsilon _{l}^{\,c\, =\, 0.6}$',...
                '$\epsilon _{l}^{\,c\, =\, 0.2}$','$k_{\epsilon}$','Location','NorthWest');
set(lgnd,'Interpreter','Latex'); xticks([0 0.25 0.5 0.75 1.0]);

% export_fig(fig,figname);
% export_fig(fig,[figname,'.pdf'],'-dpdf');

fig = figure('Color','White','Position',fig_nompossml+[fig_nompossml(3) 0 0 0]); 
figname = [fig_loc,'estd_error'];

plot(RR{2}.Qgd_est/min(RR{2}.Qgd_est),RR{2}.e_L_est); hold on;
xlabel('$Q_{g}/\dot{Q}_{g}$'); ylabel('$\epsilon _{l}, \%$');
plot([0.01 0.99],10*[1 1],'k--');
lgnd = legend('$\epsilon _{l}^{\,c\, =\, 1.0}$','$\epsilon _{l}^{\,c\, =\, 0.6}$',...
                '$\epsilon _{l}^{\,c\, =\, 0.2}$','$k_{\epsilon}$','Location','NorthWest');
set(lgnd,'Interpreter','Latex'); xticks([0 0.25 0.5 0.75 1.0]);
axis(xs); grid on;


% export_fig(fig,figname);
% export_fig(fig,[figname,'.pdf'],'-dpdf');

%%
fig = figure('Color','White','Position',[100 100 750 700]); 
figname = [fig_loc,'nominal_result_110'];

subplot(221)
plot(RR.Qgd(:,1:2:5)/min(RR.Qgd(:,1:2:5)),RR.D_Eg(:,1:2:5)); hold on;
plot(RR.Qgd(:,1:2:5)/min(RR.Qgd(:,1:2:5)),RR.D_Et(:,1:2:5),'--'); hold on;
xlabel('Qg/dot(Qg)');
ylabel('Energy/day');
xs = axis; grid on;
legend('Eg, 1.0','Eg, 0.6','Eg, 0.2','Et, 1.0','Et, 0.6','Et, 0.2','Location','NorthWest');
title('Measured Utility');

subplot(222)
plot(RR.Qgd_est(:,1:2:5)/min(RR.Qgd_est(:,1:2:5)),RR.D_Eg_est(:,1:2:5)); hold on;
plot(RR.Qgd_est(:,1:2:5)/min(RR.Qgd_est(:,1:2:5)),RR.D_Et_est(:,1:2:5),'--'); hold on;
xlabel('Qg/dot(Qg)');
ylabel('Energy/day');
axis(xs); grid on;
legend('Eg, 1.0','Eg, 0.6','Eg, 0.2','Et, 1.0','Et, 0.6','Et, 0.2','Location','NorthWest');
title('Estimated Utility');

subplot(223)
plot(RR.Qgd/min(RR.Qgd),RR.e_L); hold on;
xlabel('Qg/dot(Qg)');
ylabel('eps_L, %');
xs = axis; grid on; 
plot([0.01 0.99],10*[1 1],'k--');
legend('1.0','0.8','0.6','0.4','0.2','Location','NorthWest');
title('Measured Error');

subplot(224)
plot(RR.Qgd_est/min(RR.Qgd_est),RR.e_L_est); hold on;
xlabel('Qg/dot(Qg)');
ylabel('eps_L %');
legend('1.0','0.8','0.6','0.4','0.2','Location','NorthWest');
title('Estimated Error');
xs = axis; grid on;
plot([0.01 0.99],10*[1 1],'k--');

% export_fig(fig,figname);
% export_fig(fig,[figname,'.pdf'],'-dpdf');




