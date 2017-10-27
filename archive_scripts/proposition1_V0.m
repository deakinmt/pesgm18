% This script creates the `measures of efficiency' in Figure 3 - thermal
% efficiency, and the power factors at generator and feeder head.

clear all; close all; 
fig_loc = 'C:\Users\chri3793\Documents\DPhil\malcolm_updates\wc171023\figures\';

% addpath('pesgm_funcs');

set(0,'defaulttextinterpreter','latex');
set(0,'defaultaxesfontsize',12);
fig_nompos = [100 100 550 320];
%% First off: Maximum losses, assuming 100% at MLIMPT

fig = figure('Color','White');
fig_name = [fig_loc,'mlimpt_loss'];

Vg = 1.06;
v0_lin = linspace(0.9,1.1,1000);
theta = linspace(0.01,(pi/2) - 0.01,300);

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

[~,~,Plb,~] = V2_Sl_calc( Snb,Z,V0,'p' );
[~,~,Plp,~] = V2_Sl_calc( Snp,Z,V0,'p' );

dL = Plp - Plb;
dP = real(Snp - Snb);

dLdP = dL./dP;

cc = contourf(lz,V0,100*dLdP ); hold on; clabel(cc)

xs = axis;
lz_min = 0.14;
plot(lz_min*[1 1],xs(3:4),'k--');

set(gca,'xscale','log'); 
xlabel('$\lambda _{z} = R/X$'); ylabel('$V_{0}$');

% export_fig(fig,fig_name);
% export_fig(fig,[fig_name,'.pdf'],'-dpdf');

%% The minimum error: differential loss (numerical version)
Pnb = real(Snb);
dlta = 1e-5;
Pnbb = Pnb + dlta*dP;
[ ~,Qnbb ] = pq_pv( Pnbb,Z,Vg,V0 );
Snbb = Pnbb+ 1i*Qnbb;

[~,~,Plbb,~] = V2_Sl_calc( Snbb,Z,V0,'p' );

dLbb = Plbb - Plb;
dPnbb = real(Snbb - Snb);

dLbbdPbb = 100*dLbb./dPnbb;
cc = contourf(lz,V0,dLbbdPbb ); hold on; clabel(cc)
xs = axis;
V_min = 1.008;
lz_max = 0.63;
plot(xs(1:2),V_min*[1 1],'k-.');
plot(lz_max*[1 1],xs(3:4),'k');

set(gca,'xscale','log'); 
xlabel('$\lambda _{z}$'); ylabel('$V_{0}$');

%% Curvature
d2L = dLbbdPbb - Plb;

d2LdP = d2L./dPnbb;
cc = contourf(lz,V0,log(abs(d2LdP)) ); hold on; clabel(cc)
set(gca,'xscale','log'); 
xlabel('$\lambda _{z}$'); ylabel('$V_{0}$');



%% The minimum error: differential loss (analytic version)

dQgdPg = ( Pnb.*(aZ.^2) - (R.*(Vg.^2)) )./ ...
            sqrt( (X.*(Vg.^2)).^2 - (aZ.^2).*( ...
                    (aZ.*Pnb).^2 - 2*Pnb.*R.*(Vg.^2) + (Vg.^2).*(Vg.^2 - V0.^2)...
                                                    )  );

dPldPg = ((2*R)./(abs(Z).^2)) .*(R + (X.*dQgdPg) );
dPldPg  = dPldPg + 0./(imag(dPldPg)==0);

cc = contourf(lz,V0,dPldPg ); hold on; clabel(cc); set(gca,'xscale','log'); 

xlabel('$\lambda _{z} = R/X$'); ylabel('$V_{0}$');


%%

cc = contour(lz,nu,dLbbdPbb,'b' ); hold on; clabel(cc); hold on;
cc = contour(lz,nu,100*dLdP,'g' ); hold on; clabel(cc)
contour(lz,nu,(imag(Snb)==0),'k');
set(gca,'xscale','log'); 
xlabel('$\lambda _{z}$'); ylabel('$V_{0}/V_{+}$');
















