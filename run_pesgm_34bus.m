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

F.pg_ssc = linspace(-1e-3,0.2,100);
F.qg_ssc = linspace(-0.4,1e-3,100);
F.n = 40;

% F.pg_ssc = linspace(-1e-3,0.2,300);
% F.qg_ssc = linspace(-0.4,1e-3,500);
% F.n = 20;

% F.V0 = 1.10;
% F.Ps0_k = [1.0 0.75 0.5];

F.V0 = 1.05;
F.Ps0_k = [1.0 0.7 0.4];

% F.V0 = 1.00;
% F.Ps0_k = [1.0 0.75 0.5];

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
fig = figure('Color','White','Position',[100 100 1000 900]); 
figname = [fig_loc,'nominal_result_'];

subplot(221)
plot(RR.Qn/min(RR.Qn),RR.D_En); hold on;
plot(RR.Qn/min(RR.Qn),RR.D_Et,'--'); hold on;

plot(RR.Qn/min(RR.Qn),RR.D_En); hold on;
plot(RR.Qn/min(RR.Qn),RR.D_Et,'--'); hold on;

xlabel('Qn/dot(Qn)');
ylabel('Utility, (energy/day)');
xs = axis; grid on;

subplot(222)
plot(RR.Qn_est/min(RR.Qn_est),RR.D_En_est); hold on;
plot(RR.Qn_est/min(RR.Qn_est),RR.D_Et_est,'--'); hold on;

xlabel('Qn/dot(Qn)');
ylabel('Utility, estimated');
axis(xs); grid on;

subplot(224)
% plot(RR.Qn_est/min(RR.Qn_est),RR.e_L_est);
plot(-RR.Qn_est,RR.e_L_est);

xlabel('Qn/dot(Qn)');
ylabel('eps_L estimated, %');
grid on;
xs = axis;

subplot(223)
% plot(RR.Qn/min(RR.Qn),RR.e_L); hold on;
plot(-RR.Qn,RR.e_L); hold on;

xlabel('Qn/dot(Qn)');
ylabel('eps_L, %');
axis(xs);
xs = axis; grid on; 

% export_fig(fig,figname);
% export_fig(fig,[figname,'.pdf'],'-dpdf');




