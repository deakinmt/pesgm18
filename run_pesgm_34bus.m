close all; clear all; clc;
% cd('C:\Users\chri3793\Documents\MATLAB\DPhil\pesgm18_mtlb')
cd('C:\Users\Matt\Documents\MATLAB\DPhil\pesgm18_mtlb');
addpath('pesgm_funcs');
% fig_loc = [pwd,'\figures\'];
% fig_loc = 'C:\Users\chri3793\Documents\MATLAB\DPhil\pesgm18_mtlb\figures\';
fig_loc = 'C:\Users\Matt\Desktop\wc171023\figures\';
%%
% cd('C:\Users\chri3793\Documents\MATLAB\DPhil\pesgm18_mtlb')
cd('C:\Users\Matt\Documents\MATLAB\DPhil\pesgm18_mtlb')
% INPUTS TO DETERMINE BEHAVIOUR -------------------
F.Vp = 1.06; %pu

% F.pg_ssc = linspace(-1e-6,0.2,100);
% F.qg_ssc = linspace(-0.4,1e-6,100);
% F.n = 30;

F.pg_ssc = linspace(-1e-6,0.2,400);
F.qg_ssc = linspace(-0.4,1e-6,600);
F.n = 80;

% F.V0 = 1.00;
F.V0 = 1.10;
% F.V0 = 1.00;
F.Ps0_k = [1.0 0.8 0.6 0.4 0.2];

% remain the same:
F.SRC = 'SOURCEBUS';

F.filename = '\opendss_models\34Bus\ieee34Mod1_fxd';
F.feeder = '34bus'; % for calc_ww and T_ratios
F.tr_buses = {{'888','832'};{'832','852'};{'850','814'};{'800','SOURCEBUS'}};
F.fig_nompos = [200 150 450 340];

SCD = load('solar_curves_dataset.mat','Ps');
F.Ps0 = SCD.Ps; %pu
% F.Ps = 2.7*SCD.Ps; %pu
F.Qmaxpu = 4;

% First results figure: 834 bus feeder voltage
figname1 = [fig_loc,'34bus_pgp0'];
F.NUT = '834';

pl_options={'X'};
[ RR,figs ] = run_pesgm_feeder( F,pl_options );

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




