% This script creates the `measures of efficiency' in Figure 3 - thermal
% efficiency, and the power factors at generator and feeder head.

clear all; close all; 
fig_loc = 'C:\Users\chri3793\Documents\DPhil\malcolm_updates\wc171016\figures\';

set(0,'defaulttextinterpreter','latex');
set(0,'defaultaxesfontsize',12);
fig_nompos = [100 100 550 320];
%% 

Vg = 1.0;
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

%%
fig = figure('Color','White','Position',fig_nompos+[0 0 0 300]);
fig_name=[fig_loc,'results1'];

dPgdQg = -((lv.^2).*lx)./sqrt( (lv.^2) - ( (lv.^4).*(lx.^2) ) );
dPtdQg = (lx.^2).*( -2 + ( lz + lz.^-1 ).*dPgdQg );

dPtdQg = dPtdQg.*(imag(dPtdQg)==0);
dPgdQg = dPgdQg.*(imag(dPgdQg)==0);

subplot(221)
cc = contourf(lz,lv,log10(dPtdQg*-1) );
set(gca,'xscale','log'); clabel(cc)
legend('log10( -dPtdQn(Qn = 0) )','Location','NorthWest');
xlabel('$\lambda _{z}$'); ylabel('$\lambda _{v}$');

subplot(222)
cc = contourf(lz,lv,log10((dPgdQg)*-1) );
set(gca,'xscale','log'); clabel(cc)
legend('log10( -dPtdQn(Qn = 0) )','Location','NorthWest');
xlabel('$\lambda _{z}$'); ylabel('$\lambda _{v}$');

subplot(223)
cc = contourf(lz,lv,log10(abs(dPgdQg./dPtdQg)),(-2:0.1:2));
set(gca,'xscale','log'); clabel(cc);
legend('log10( -dPtdQn(Qn = 0) )','Location','NorthWest');
xlabel('$\lambda _{z}$'); ylabel('$\lambda _{v}$');

subplot(224)
cc = contourf(lz,lv,log10( abs((dPtdQg-dPgdQg)./dPtdQg) ),(-1:0.25:3) );
set(gca,'xscale','log'); clabel(cc);
legend('log10( -dPtdQn(Qn = 0) )','Location','NorthWest');
xlabel('$\lambda _{z}$'); ylabel('$\lambda _{v}$');


% subplot(224)
% cc = contourf(lz,lv,sign(dPgdQg-dPtdQg) );
% set(gca,'xscale','log'); clabel(cc);
% legend('log10( -dPtdQn(Qn = 0) )','Location','NorthWest');
% xlabel('$\lambda _{z}$'); ylabel('$\lambda _{v}$');

% grid on;

% export_fig(fig,fig_name);
% export_fig(fig,[fig_name,'.pdf'],'-dpdf');


%%
[ Pnb,~ ] = pq_pv( 0,Z,Vg,V0 );
Pnb = Pnb.*(imag(Pnb)==0);
[V2,Sl,Pl,Ql] = V2_Sl_calc( Pnb,Z,V0,'p' );

figure;
cc = contourf(lz,lv,Pl,(0.0001:0.005:0.05));
set(gca,'xscale','log'); clabel(cc);


%%

dQgdPg = -((lv.^2).*lr)./sqrt( (lv.^2) - ( (lv.^4).*(lr.^2) ) );


dPldPg = ((2*R)./(abs(Z).^2)).*(R + (X.*dQgdPg) );
dPldPg  = dPldPg + 0./(imag(dPldPg)==0);

dP0dPg = 1 - dPldPg;

% subplot(121)
% cc = contourf(lz,V0,log10(abs(dP0dPg)),(-1:0.1:1) );
cc = contourf(lz,V0,real(dP0dPg),(0:0.1:1.5));
% cc = contourf(lz,V0,log10(abs(dPldPg)) );
set(gca,'xscale','log'); clabel(cc)
legend('log10( -dPtdQn(Qn = 0) )','Location','NorthWest');
xlabel('$\lambda _{z}$'); ylabel('$V_{0}$');

grid on;

% subplot(122)
% cc = contourf(lz,V0,sign(dPldPg) );
% set(gca,'xscale','log'); clabel(cc)
% legend('log10( -dPtdQn(Qn = 0) )','Location','NorthWest');
% xlabel('$\lambda _{z}$'); ylabel('$V_{0}$');

















