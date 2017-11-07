clear all; close all; 
fig_loc = 'C:\Users\chri3793\Documents\DPhil\pesgm18\pesgm18_paper\figures\';
addpath('pesgm_funcs');

set(0,'defaulttextinterpreter','latex');
set(0,'defaultaxesfontsize',12);
fig_nompos = [100 100 550 300];

Vg = 1.06;
v0_lin = linspace(0.9,1.1,300);
theta = linspace(0.01,(pi/2) - 0.01,300);
% v0_lin = linspace(0.9,1.1,3000); % higher resolution if needs be
% theta = linspace(0.01,(pi/2) - 0.01,1000); % higher resolution if needs be

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
%% Plot k'
fig = figure('Color','White','Position',fig_nompos);
fig_name = [fig_loc,'k_dot'];

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

cc = contourf(lz,V0,eps_L,[0 10 (20:20:300)]); hold on; 
set(gca,'xscale','log');
clabel(cc)
% clabel(cc,'manual'); % to make it prettier

xlbl = xlabel('$R/X$');  ylbl = ylabel('\textit{V$_{t}$}');
set(xlbl,'FontSize',14); set(ylbl,'FontSize',14);
colormap parula

lgnd = legend('$k''$, \%','Location','SouthWest');
set(lgnd,'Interpreter','Latex','FontSize',14,'Box','Off');
% export_fig(fig,fig_name); % NB the PDF version doesn't work!
%% Plot knom
fig = figure('Color','White','Position',fig_nompos);
fig_name = [fig_loc,'k_bar'];

Kpq = (-X*(Vg.^2))./sqrt( ((R.*(Vg.^2)).^2) - (((aZ.*Vg).^2).*(Vg.^2 - V0.^2)) );
Knum = 2.*R.*(R.*Kpq + X);
Kb = Knum./( Knum - ((aZ.^2).*Kpq) );

Kb  = Kb + 0./(imag(Kb)==0);

cc = contourf(lz,V0,-real(100*Kb),100*(-0.2:0.05:1.15)); hold on; 
set(gca,'xscale','log'); 
clabel(cc);
% clabel(cc,'manual');
xlbl = xlabel('$R/X$');  ylbl = ylabel('\textit{V$_{t}$}');
set(xlbl,'FontSize',14); set(ylbl,'FontSize',14);

lgnd = legend('$k_{nom}\,$ , \%','Location','SouthWest');
set(lgnd,'Interpreter','Latex','FontSize',14,'Box','Off');
colormap parula
% export_fig(fig,fig_name);





