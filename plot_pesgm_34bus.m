close all; clear all;

set(0,'defaulttextinterpreter','latex');
set(0,'defaultaxesfontsize',14);
fig_nompossml = [100 100 300 350];
fig_nompos = [100 100 420 350];

data_name = 'run_pesgm_34bus_rslt.mat'
load(data_name);

% fig_loc = [pwd,'\figures\'];
% fig_loc = 'C:\Users\chri3793\Documents\DPhil\pesgm18\pesgm18_paper\figures\';
fig_loc = 'C:\Users\Matt\Documents\DPhil\pesgm18\pesgm18_paper\figures\';



%% UPE
% fig = figure('Color','White','Position',fig_nompos); % uncomment for larger
fig = figure('Color','White','Position',fig_nompossml);
figname = [fig_loc,'UPE'];

plot(RR{1}.DPg/(max(RR{1}.DPg)),RR{1}.eps_P*100); hold on;
plot(RR{2}.DPg/(max(RR{2}.DPg)),RR{2}.eps_P*100);
plot(RR{3}.DPg/(max(RR{3}.DPg)),RR{3}.eps_P*100);

plot(RR{1}.DPg_estd/max(RR{1}.DPg_estd),RR{1}.eps_P_estd*100,'--');
plot(RR{2}.DPg_estd/max(RR{2}.DPg_estd),RR{2}.eps_P_estd*100,'--');
plot(RR{3}.DPg_estd/max(RR{3}.DPg_estd),RR{3}.eps_P_estd*100,'--');

lgnd = legend('$\epsilon _{P,\, msrd}^{Vt=1.04}$','$\epsilon _{P,\, msrd}^{Vt=1.08}$',...
               '$\epsilon _{P,\, msrd}^{Vt=1.12}$','$\epsilon _{P,\, estd}^{Vt=1.04}$',...
               '$\epsilon _{P,\, estd}^{Vt=1.08}$','$\epsilon _{P,\, estd}^{Vt=1.12}$',...
                                                    'Location','NorthWest');
% set(lgnd,'Interpreter','Latex','Fontsize',15); % uncomment for larger
set(lgnd,'Interpreter','Latex');
xlabel('$\Delta P_{g}/\Delta P_{g}^{\prime}$');
ylabel('$\epsilon_{P}$, \%'); grid on;
xticks([0 0.25 0.5 0.75 1.0]);
axis([0 1 -20 220]);

% export_fig(fig,figname);
% export_fig(fig,[figname,'.pdf'],'-dpdf');
%% PLOT #1 + #2: Measured and Estimated Utility versus c VERSION 2!!!!!!!
fig = figure('Color','White','Position',fig_nompossml); 
figname = [fig_loc,'eg_utility'];

idx = 4;

sb_mw = RR{idx}.sbase*1e-3;

plot( RR{idx}.Qgd/min(RR{idx}.Qgd) , RR{idx}.D_Eg*sb_mw ); hold on;
plot( RR{idx}.Qgd_est/min(RR{idx}.Qgd_est) , RR{idx}.D_Eg_est*sb_mw,'--' );
xlabel('$Q_{g}/\tilde{Q}_{g}$'); ylabel('$\Delta E_{g}$, MWh');
lgnd = legend('$\Delta E_{g,\,msrd}^{\,c\, =\, 1.0}$','$\Delta E_{g,\,msrd}^{\,c\, =\, 0.6}$',...
                '$\Delta E_{g,\,msrd}^{\,c\, =\, 0.2}$','$\Delta E_{g,\,estd}^{\,c\, =\, 1.0}$',...
                '$\Delta E_{g,\,estd}^{\,c\, =\, 0.6}$','$\Delta E_{g,\,estd}^{\,c\, =\, 0.2}$','Location','NorthWest');
set(lgnd,'Interpreter','Latex');
xticks([0 0.25 0.5 0.75 1.0]);
axis([0 1 0 46]); grid on;

% export_fig(fig,figname);
% export_fig(fig,[figname,'.pdf'],'-dpdf');
%-------------------------
fig = figure('Color','White','Position',fig_nompossml+[fig_nompossml(3) 0 0 0]); 
figname = [fig_loc,'et_utility'];

plot( RR{idx}.Qgd/min(RR{idx}.Qgd) , RR{idx}.D_Et*sb_mw); hold on;
plot( RR{idx}.Qgd_est/min(RR{idx}.Qgd_est) , RR{idx}.D_Et_est*sb_mw ,'--');
xlabel('$Q_{g}/\tilde{Q}_{g}$'); ylabel('$\Delta E_{t}$, MWh');
lgnd = legend('$\Delta E_{t,\, msrd}^{\,c\, =\, 1.0}$','$\Delta E_{t,\, msrd}^{\,c\, =\, 0.6}$',...
                '$\Delta E_{t,\, msrd}^{\,c\, =\, 0.2}$','$\Delta E_{t,\,estd}^{\,c\, =\, 1.0}$',...
                '$\Delta E_{t,\,estd}^{\,c\, =\, 0.6}$','$\Delta E_{t,\,estd}^{\,c\, =\, 0.2}$','Location','NorthWest');
set(lgnd,'Interpreter','Latex'); xticks([0 0.25 0.5 0.75 1.0]);
axis([0 1 0 46]); grid on;

% export_fig(fig,figname);
% export_fig(fig,[figname,'.pdf'],'-dpdf');

%% UEE
fig = figure('Color','White','Position',fig_nompos); 
figname = [fig_loc,'UEE'];

plot(RR{idx}.Qgd/min(RR{idx}.Qgd),RR{idx}.e_L); hold on;
plot(RR{idx}.Qgd_est/min(RR{idx}.Qgd_est),RR{idx}.e_L_est,'--');

xlabel('$Q_{g}/\tilde{Q}_{g}$'); ylabel('$\epsilon _{E}, \%$');
axis([0 1 0 140]);
xs = axis; grid on;
lgnd = legend('$\epsilon _{E,\,msrd}^{\,c\, =\, 1.0}$','$\epsilon _{E,\,msrd}^{\,c\, =\, 0.6}$',...
              '$\epsilon _{E,\,msrd}^{\,c\, =\, 0.2}$','$\epsilon _{E,\,estd}^{\,c\, =\, 1.0}$',...
              '$\epsilon _{E,\,estd}^{\,c\, =\, 0.6}$','$\epsilon _{E,\,estd}^{\,c\, =\, 0.2}$',...
                                                            'Location','NorthWest');
set(lgnd,'Interpreter','Latex','Fontsize',15); xticks([0 0.25 0.5 0.75 1.0]);
yticks([0 25 50 75 100 125]);

% export_fig(fig,figname);
% export_fig(fig,[figname,'.pdf'],'-dpdf');
%------------------

