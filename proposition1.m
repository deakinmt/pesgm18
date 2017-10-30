% This script creates the `measures of efficiency' in Figure 3 - thermal
% efficiency, and the power factors at generator and feeder head.

clear all; close all; 
fig_loc = 'C:\Users\chri3793\Documents\DPhil\pesgm18\pesgm18_paper\figures\';

cd('C:\Users\chri3793\Documents\MATLAB\DPhil\pesgm18_mtlb');addpath('pesgm_funcs');

set(0,'defaulttextinterpreter','latex');
set(0,'defaultaxesfontsize',12);
fig_nompos = [100 100 550 320];
%% First off: Maximum losses, assuming 100% at MLIMPT
fig = figure('Color','White','Position',fig_nompos);
fig_name = [fig_loc,'mlimpt_loss_2013b'];

Vg = 1.06;
v0_lin = linspace(0.9,1.1,1000);
theta = linspace(0.01,(pi/2) - 0.01,300);
% v0_lin = linspace(0.9,1.1,3000);
% theta = linspace(0.01,(pi/2) - 0.01,1000);

lz_lin = cot(theta); % lambda
lv_lin = Vg./v0_lin;
nu_lin = 1./lv_lin;

[lz,nu] = meshgrid(lz_lin,nu_lin);
lv = 1./nu;
V0 = Vg./lv;

lr = lz./sqrt(1 + (lz.^2));
lx = 1./sqrt(1 + (lz.^2));
R = lr; X = lx;
Z = R + 1i*X;
aZ = abs(Z);

[ Snb,~,Snp ] = calc_xpts( 0,0,Z,Vg,V0 );
Snb = Snb + 0./(imag(Snb)==0);

Pnp = real(Snp);
Pnb = real(Snb);

[~,~,Plb,~] = V2_Sl_calc( Snb,Z,V0,'p' );
[~,~,Plp,~] = V2_Sl_calc( Snp,Z,V0,'p' );

Ptp = (Pnp - Plp);
Ptb = (Pnb - Plb);
DPt = Ptp - Ptb;
DPn = Pnp - Pnb;

eps_L = 100*(DPn - DPt)./DPt;

% dL = Plp - Plb;
% dP = real(Snp - Snb);
% dLdP = dL./dP;
% cc = contourf(lz,V0,100*dLdP ); hold on; clabel(cc)
cc = contourf(lz,V0,eps_L,[0 10 (20:20:300)]); hold on; clabel(cc)

xs = axis;
lz_min = 0.09;
plot(lz_min*[1 1],xs(3:4),'k--');

set(gca,'xscale','log'); 
xlabel('$\lambda _{z} = R/X$'); ylabel('$V_{0}$');
colormap parula
% export_fig(fig,fig_name);
% export_fig(fig,[fig_name,'.pdf'],'-dpdf');

%% The minimum error: differential loss (analytic version)
fig = figure('Color','White','Position',fig_nompos);
fig_name = [fig_loc,'dPldPg_2013b'];

Pnb = real(Snb);

BX = (X.*(Vg.^2)).^2 - (aZ.^2).*( ...
                    (aZ.*Pnb).^2 - 2*Pnb.*R.*(Vg.^2) + (Vg.^2).*(Vg.^2 - V0.^2)...
                                                    );

dQgdPg = ( Pnb.*(aZ.^2) - (R.*(Vg.^2)) )./sqrt(BX);

dPldPg = ((2*R)./(abs(Z).^2)) .*(R + (X.*dQgdPg) );
dPldPg  = dPldPg + 0./(imag(dPldPg)==0);

[cc,hh] = contourf(lz,V0,100*dPldPg,100*(-0.15:0.05:0.7) ); hold on; 
QQQQ = 100*[-0.05,0,0.05,(0.1:0.1:0.5)];
clabel(cc,QQQQ); set(gca,'xscale','log'); 
xlabel('$\lambda _{z} = R/X$'); ylabel('$V_{0}$');
xs = axis;
V_min = 1.007;
lz_max = 0.63;
plot(xs(1:2),V_min*[1 1],'k-.');
plot(lz_max*[1 1],xs(3:4),'k');

colormap parula
% export_fig(fig,fig_name);
% export_fig(fig,[fig_name,'.pdf'],'-dpdf');