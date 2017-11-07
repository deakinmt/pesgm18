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

PF = 1 - linspace(0,2*dPF,5);

Sload = aS*exp(1i*acos(PF));

fig = figure('Color','White','Position',fig_nompossml);
figname = [fig_loc,'nonlin_S0_PF'];

for i = 1:numel(PF)
    [ ~,~,~,~,eps_P,DPg ] = estm_DE_EL( 1,1,Sload(i),Z,Vp,V0,n );
    plot(DPg./max(DPg),100*eps_P,'--'); hold on;
end

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