% This script creates the `measures of efficiency' in Figure 3 - thermal
% efficiency, and the power factors at generator and feeder head.

clear all; close all; 
fig_loc = 'C:\Users\chri3793\Documents\DPhil\malcolm_updates\wc171016\figures\';

set(0,'defaulttextinterpreter','latex');
set(0,'defaultaxesfontsize',12);
fig_nompos = [100 100 550 320];
%% 
fig = figure('Color','White','Position',fig_nompos);
fig_name=[fig_loc,'results1'];

Vg = 1.06;
v0_lin = linspace(0.9,1.2,300);
theta = linspace(0.003,(pi/2) - 0.0006,1000);

lz_lin = cot(theta); % lambda
lv_lin = Vg./v0_lin;
nu_lin = 1./lv_lin;

[lz,nu] = meshgrid(lz_lin,nu_lin);
lv = 1./nu;
V0 = Vg./lv;

lr = lz./sqrt(1 + (lz.^2));
lx = 1./sqrt(1 + (lz.^2));
Z = lr + 1i*lx;

dPgdQg = -((lv.^2).*lx)./sqrt( (lv.^2) - ( (lv.^4).*(lx.^2) ) );

dPtdQg = (lx.^2).*( -2 + ( lz + lz.^-1 ).*dPgdQg );
dPtdQg = dPtdQg.*(imag(dPtdQg)==0);

cc = contourf(lz,lv,log10(dPtdQg*-1) );
% cc = contour(lz,lv,dPtdQg*-1 );
set(gca,'xscale','log');
clabel(cc)
legend('log10( -dPtdQn(Qn = 0) )','Location','NorthWest');

xlabel('$\lambda _{z}$'); ylabel('$\lambda _{v}$');
% grid on;

export_fig(fig,fig_name);
export_fig(fig,[fig_name,'.pdf'],'-dpdf');










