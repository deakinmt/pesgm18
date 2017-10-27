close all; clear all; clc;
fig_loc = 'C:\Users\chri3793\Documents\DPhil\malcolm_updates\wc171023\figures\';

%% first plot the magnitude of the curvature at Pg

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

fig = figure('Color','White');
fig_name = [fig_loc,'curvature_plots_d2LdP2'];

Pnb = real(Snb);

BX = (X.*(Vg.^2)).^2 - (aZ.^2).*( ...
                    (aZ.*Pnb).^2 - 2*Pnb.*R.*(Vg.^2) + (Vg.^2).*(Vg.^2 - V0.^2)...
                                                    );
d2Pl_dPg2 = ((aZ.^2)./(BX.^1.5)).*( BX + ( Pnb.*(aZ.^2) - R*(Vg.^2) ).^2 );

[cc,hh] = contourf(lz,V0,log10(d2Pl_dPg2),(-1:0.5:8) ); hold on; 
clabel(cc); set(gca,'xscale','log'); 
xlabel('$\lambda _{z} = R/X$'); ylabel('$V_{0}$');
title('Log10(d2L/dP2)');

% export_fig(fig,fig_name);
% export_fig(fig,[fig_name,'.pdf'],'-dpdf');

%% PLOT curvature versus Pg
% Zexpn0 = pi*[1/2.00001 1/2.1 1/3 1/4 1/6 1/20 1e-10];
% Zexpn0 = pi*[1/2.00001 1/3 1e-10];

Zexpn0 = pi*[1/2.1 1/2.5 1/4 1e-2];
V00 = 1.00;
Vg0 = 1.06;

Pg0 = (-1.3:0.001:3);
Z0 = exp(1i*Zexpn0);
lz = real(Z0)./imag(Z0)

fig = figure('Color','White','Position',[100 100 950 650]);
PP = [];
for i = 1:numel(Z0)
    R0 = real(Z0(i));
    X0 = imag(Z0(i));
    aZ0 = abs(Z0(i));
    
    [ Snb0,~,Snp0 ] = calc_xpts( 0,0,Z0(i),Vg0,V00 )
    

    BX0 = (X0.*(Vg0.^2)).^2 - (aZ0.^2).*( ...
                        (aZ0.*Pg0).^2 - 2*Pg0.*R0.*(Vg0.^2) + (Vg0.^2).*(Vg0.^2 - V00.^2)...
                                                        );

    d2Pl_dPg20 = ((aZ0.^2)./(BX0.^1.5)).*( BX0 + ( Pg0.*(aZ0.^2) - R0*(Vg0.^2) ).^2 );
    d2Pl_dPg20 = d2Pl_dPg20 + 0./(imag(d2Pl_dPg20)==0)*(1 + 1i);
    PP = [PP;plot(Pg0,d2Pl_dPg20)]; grid on; hold on;
    xs = axis;
    Snb0 = Snb0 + 0./(imag(Snb0)==0);
    
    plot(real(Snb0)*[1 1],xs(3:4),'r');
    plot(real(Snp0)*[1 1],xs(3:4),'k--');
end
axis([-1.2 2.2 0.25 10]);
xlabel('Pn');
ylabel('d2L/dP2');

legend(PP,'R/X = 0.075','R/X = 0.32','R/X = 1','R/X = 30')
title('V0 = 1.00, Vg = 1.06');
fig_name = [fig_loc,'curvature_plots100106'];

export_fig(fig,fig_name);
export_fig(fig,[fig_name,'.pdf'],'-dpdf');


