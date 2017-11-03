close all; clear all; clc;
% cd('C:\Users\chri3793\Documents\MATLAB\DPhil\pesgm18_mtlb')
% cd('C:\Users\Matt\Documents\MATLAB\DPhil\pesgm18_mtlb');
addpath('pesgm_funcs');
% fig_loc = [pwd,'\figures\'];
% fig_loc = 'C:\Users\chri3793\Documents\MATLAB\DPhil\pesgm18_mtlb\figures\';
fig_loc = 'C:\Users\Matt\Desktop\wc171023\figures\';

set(0,'defaulttextinterpreter','latex');
set(0,'defaultaxesfontsize',14);
fig_nompossml = [100 100 300 350];
%%
% INPUTS TO DETERMINE BEHAVIOUR -------------------
% C = {[1.0 0.6 0.2],0.4,0.4,0.4,0.4};
% VT = {1.05,1.10,1.06,1.02,0.98};

C = {1.0,1.0,1.0};
VT = {1.04,1.08,1.12};

% WARNING! takes 3-4 minutes per run.
% F.pg_ssc = linspace(-1e-6,0.2,400);
% F.qg_ssc = linspace(-0.4,1e-6,600);
F.n = 80;
F.pg_ssc = linspace(-1e-6,0.2,100);
F.qg_ssc = linspace(-0.4,1e-6,100);
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


%%
fig = figure('Color','White','Position',fig_nompossml); 
plot(RR{1}.Qln/(min(RR{1}.Qln)),RR{1}.eps_P*100); hold on;
plot(RR{2}.Qln/(min(RR{2}.Qln)),RR{2}.eps_P*100);
plot(RR{3}.Qln/(min(RR{3}.Qln)),RR{3}.eps_P*100); hold on;
lgnd = legend('$\epsilon _{P}^{Vt=1.04}$','$\epsilon _{P}^{Vt=1.08}$',...
                        '$\epsilon _{P}^{Vt=1.12}$','Location','NorthWest');
set(lgnd,'Interpreter','Latex');                    
xlabel('$Q_{g}/Q_{g}^{\prime}$');
ylabel('$\epsilon_{P}$, \%'); grid on;
xticks([0 0.25 0.5 0.75 1.0]);
% export_fig(fig,figname);
% export_fig(fig,[figname,'.pdf'],'-dpdf');

fig = figure('Color','White','Position',fig_nompossml+[fig_nompossml(3) 0 0 0]); 
plot(RR{1}.Qgd_est/min(RR{1}.Qgd_est),RR{1}.eps_P_estd*100); hold on;
plot(RR{2}.Qgd_est/min(RR{2}.Qgd_est),RR{2}.eps_P_estd*100);
plot(RR{3}.Qgd_est/min(RR{3}.Qgd_est),RR{3}.eps_P_estd*100);
lgnd = legend('$\epsilon _{P}^{Vt=1.04}$','$\epsilon _{P}^{Vt=1.08}$',...
                        '$\epsilon _{P}^{Vt=1.12}$','Location','NorthWest');
set(lgnd,'Interpreter','Latex');
xlabel('$Q_{g}/Q_{g}^{\prime}$');
ylabel('$\epsilon_{P}$, \%'); grid on;
xticks([0 0.25 0.5 0.75 1.0]);
% export_fig(fig,figname);
% export_fig(fig,[figname,'.pdf'],'-dpdf');
%% PLOT #1 + #2: Measured and Estimated Utility versus c
fig = figure('Color','White','Position',fig_nompossml); 
figname = [fig_loc,'msrd_c_utility'];

sb_mw = RR{1}.sbase*1e-3;

plot( RR{1}.Qgd/min(RR{1}.Qgd) , RR{1}.D_Eg*sb_mw ); hold on;
plot( RR{1}.Qgd/min(RR{1}.Qgd) , RR{1}.D_Et*sb_mw ,'--'); hold on;
xlabel('$Q_{g}/\dot{Q}_{g}$'); ylabel('$\Delta E_{(\cdot)}$, MWh');
lgnd = legend('$\Delta E_{g}^{\,c\, =\, 1.0}$','$\Delta E_{g}^{\,c\, =\, 0.6}$',...
                '$\Delta E_{g}^{\,c\, =\, 0.2}$','$\Delta E_{t}^{\,c\, =\, 1.0}$',...
                '$\Delta E_{t}^{\,c\, =\, 0.6}$','$\Delta E_{t}^{\,c\, =\, 0.2}$','Location','NorthWest');
set(lgnd,'Interpreter','Latex');
xticks([0 0.25 0.5 0.75 1.0]);
axis([0 1 0 46]);
xs = axis; grid on;

% export_fig(fig,figname);
% export_fig(fig,[figname,'.pdf'],'-dpdf');

fig = figure('Color','White','Position',fig_nompossml+[fig_nompossml(3) 0 0 0]); 
figname = [fig_loc,'estd_c_utility'];

plot( RR{1}.Qgd_est/min(RR{1}.Qgd_est) , RR{1}.D_Eg_est*sb_mw ); hold on;
plot( RR{1}.Qgd_est/min(RR{1}.Qgd_est) , RR{1}.D_Et_est*sb_mw ,'--'); hold on;
xlabel('$Q_{g}/\dot{Q}_{g}$'); ylabel('$\Delta E_{(\cdot)}$, MWh');
lgnd = legend('$\Delta E_{g}^{\,c\, =\, 1.0}$','$\Delta E_{g}^{\,c\, =\, 0.6}$',...
                '$\Delta E_{g}^{\,c\, =\, 0.2}$','$\Delta E_{t}^{\,c\, =\, 1.0}$',...
                '$\Delta E_{t}^{\,c\, =\, 0.6}$','$\Delta E_{t}^{\,c\, =\, 0.2}$','Location','NorthWest');
set(lgnd,'Interpreter','Latex'); xticks([0 0.25 0.5 0.75 1.0]);
axis(xs); grid on;

% export_fig(fig,figname);
% export_fig(fig,[figname,'.pdf'],'-dpdf');

%% PLOT #3 + #4: Measured and estimated Error versus c
fig = figure('Color','White','Position',fig_nompossml); 
figname = [fig_loc,'msrd_c_error'];

plot(RR{1}.Qgd/min(RR{1}.Qgd),RR{1}.e_L); hold on;
xlabel('$Q_{g}/\dot{Q}_{g}$'); ylabel('$\epsilon _{l}, \%$');
axis([0 1 0 140]);
xs = axis; grid on;
plot([0.01 0.99],10*[1 1],'k--');
lgnd = legend('$\epsilon _{l}^{\,c\, =\, 1.0}$','$\epsilon _{l}^{\,c\, =\, 0.6}$',...
                '$\epsilon _{l}^{\,c\, =\, 0.2}$','$k_{\epsilon}$','Location','NorthWest');
set(lgnd,'Interpreter','Latex'); xticks([0 0.25 0.5 0.75 1.0]);
yticks([0 25 50 75 100 125]);

% export_fig(fig,figname);
% export_fig(fig,[figname,'.pdf'],'-dpdf');

fig = figure('Color','White','Position',fig_nompossml+[fig_nompossml(3) 0 0 0]); 
figname = [fig_loc,'estd_c_error'];

plot(RR{1}.Qgd_est/min(RR{1}.Qgd_est),RR{1}.e_L_est); hold on;
xlabel('$Q_{g}/\dot{Q}_{g}$'); ylabel('$\epsilon _{l}, \%$');
plot([0.01 0.99],10*[1 1],'k--');
lgnd = legend('$\epsilon _{l}^{\,c\, =\, 1.0}$','$\epsilon _{l}^{\,c\, =\, 0.6}$',...
                '$\epsilon _{l}^{\,c\, =\, 0.2}$','$k_{\epsilon}$','Location','NorthWest');
set(lgnd,'Interpreter','Latex'); xticks([0 0.25 0.5 0.75 1.0]);
axis(xs); grid on;
yticks([0 25 50 75 100 125]);

% export_fig(fig,figname);
% export_fig(fig,[figname,'.pdf'],'-dpdf');

%% PLOT 5 + 6: Measured and estimated error versus Vt
fig = figure('Color','White','Position',fig_nompossml); 
figname = [fig_loc,'msrd_vt_error'];

for i = 2:numel(RR)
    plot(RR{i}.Qgd/min(RR{i}.Qgd),RR{i}.e_L); hold on;
end    
xlabel('$Q_{g}/\dot{Q}_{g}$'); ylabel('$\epsilon _{l}, \%$');
axis([0 1 -7 120]);
xs = axis; grid on;
plot([0.01 0.99],10*[1 1],'k--');
lgnd = legend('$\epsilon _{l}^{\,Vt\, =\, 1.10}$','$\epsilon _{l}^{\,Vt\, =\, 1.06}$',...
                '$\epsilon _{l}^{\,Vt\, =\, 1.02}$','$\epsilon _{l}^{\,Vt\, =\, 0.98}$',...
                '$k_{\epsilon}$','Location','NorthWest');
set(lgnd,'Interpreter','Latex'); xticks([0 0.25 0.5 0.75 1.0]);

% export_fig(fig,figname);
% export_fig(fig,[figname,'.pdf'],'-dpdf');

fig = figure('Color','White','Position',fig_nompossml+[fig_nompossml(3) 0 0 0]); 
figname = [fig_loc,'estd_vt_error'];

for i = 2:numel(RR)
    plot(RR{i}.Qgd_est/min(RR{i}.Qgd_est),RR{i}.e_L_est); hold on;
end    
xlabel('$Q_{g}/\dot{Q}_{g}$'); ylabel('$\epsilon _{l}, \%$');
axis(xs); grid on;
plot([0.01 0.99],10*[1 1],'k--');
lgnd = legend('$\epsilon _{l}^{\,Vt\, =\, 1.10}$','$\epsilon _{l}^{\,Vt\, =\, 1.06}$',...
                '$\epsilon _{l}^{\,Vt\, =\, 1.02}$','$\epsilon _{l}^{\,Vt\, =\, 0.98}$',...
                '$k_{\epsilon}$','Location','NorthWest');
set(lgnd,'Interpreter','Latex'); xticks([0 0.25 0.5 0.75 1.0]);

% export_fig(fig,figname);
% export_fig(fig,[figname,'.pdf'],'-dpdf');



