close all; clear all; clc;

addpath('pesgm_funcs');
fig_loc = 'C:\Users\Matt\Documents\DPhil\pesgm18\pesgm18_paper\figures\';
set(0,'defaulttextinterpreter','latex');
set(0,'defaultaxesfontsize',14);
fig_nompossml = [100 100 300 350];

aZ = 0.203; 
lz = 1.85;
Z = aZ*exp(1i*acot(lz));
Vp = 1.06;
V0 = 1.04;
n = 500;
%%
th_ld = 9.4*pi/180;
dPF = 1 - cos(th_ld);
aS = 0.72;

PF0 = 1;
PF1 = 1 - 0.5*dPF
PF2 = 1 - dPF
PF3 = 1 - 1.5*dPF
PF4 = 1 - 2*dPF

Sload0 = aS*exp(1i*acos(PF0));
Sload1 = aS*exp(1i*acos(PF1));
Sload2 = aS*exp(1i*acos(PF2));
Sload3 = aS*exp(1i*acos(PF3));
Sload4 = aS*exp(1i*acos(PF4));

% Sload = 0.72*exp(1i*acos(PF0));
% Sload_new = 0.72*exp(1i*acos(PF_new));

[ ~,~,~,~,eps_P0,DPg0 ] = estm_DE_EL( 1,1,Sload0,Z,Vp,V0,n );
[ ~,~,~,~,eps_P1,DPg1 ] = estm_DE_EL( 1,1,Sload1,Z,Vp,V0,n );
[ ~,~,~,~,eps_P2,DPg2 ] = estm_DE_EL( 1,1,Sload2,Z,Vp,V0,n );
[ ~,~,~,~,eps_P3,DPg3 ] = estm_DE_EL( 1,1,Sload3,Z,Vp,V0,n );
[ ~,~,~,~,eps_P4,DPg4 ] = estm_DE_EL( 1,1,Sload4,Z,Vp,V0,n );

fig = figure('Color','White','Position',fig_nompossml);
figname = [fig_loc,'nonlin_S0_PF'];
plot(DPg0./max(DPg0),100*eps_P0,'--'); hold on;
plot(DPg1./max(DPg1),100*eps_P1,'--');
plot(DPg2./max(DPg2),100*eps_P2,'--');
plot(DPg3./max(DPg3),100*eps_P3,'--');
plot(DPg4./max(DPg4),100*eps_P4,'--');

% plot(DPg1./max(DPg1),100*eps_P1); hold on;
% plot(DPg_er./max(DPg_er),100*eps_P_er);

xlabel('$\Delta P_{g}/\Delta P_{g}^{\prime}$');
ylabel('$\epsilon_{P}$, \%');
xticks([0 0.25 0.5 0.75 1.0]);
axis([0 1 -20 220]); %same as UPE plot

lgnd = legend('$\epsilon_{P,\,estd}^{PF=1.00}$','$\epsilon_{P,\,estd}^{PF=0.993}$',...
              '$\epsilon_{P,\,estd}^{PF=0.987}$','$\epsilon_{P,\,estd}^{PF=0.980}$',...
              '$\epsilon_{P,\,estd}^{PF=0.973}$','Location','NorthWest');
set(lgnd,'Interpreter','Latex');
grid on;

% export_fig(fig,figname);
% export_fig(fig,[figname,'.pdf'],'-dpdf');