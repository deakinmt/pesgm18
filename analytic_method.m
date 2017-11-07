clear all; close all; clc

addpath('pesgm_funcs');

fig_loc = 'C:\Users\Matt\Documents\DPhil\pesgm18\pesgm18_paper\figures\';

set(0,'defaulttextinterpreter','latex');
set(0,'defaultaxesfontsize',14);
fig_nompos = [100 100 550 320];
%% 
Vp = 1.00;
V0 = 0.9;
Z = 0.203*exp(1i*acot(1.85)); %see PSCC paper
S0 = 0.72*exp(1i*9.4*pi/180);
% S0 = 0;
Qgen = 0.7;

ds = 0.003;

pg = 5*(-0.2:ds:0.75);
qg = 5*(-0.6:ds:0.15);

[Pg,Qg] = meshgrid(pg,qg);

Sg = Pg + 1i*Qg;

[~,Sl,~,~] = V2_Sl_calc( Sg,Z,V0,'p' );
[ Se,Sg,Snb,Snd,Snp ] = calc_s1( Sg,S0,Qgen(1),Z,Vp,V0 );
Pgh = 0.8*real(Snp);


%------------------ PLOT
fig_name = [fig_loc,'analytic_method'];
fig = figure('Color','White');
contour(Pg,Qg,Pg-real(Sl)); axis equal; hold on;
[ ~,qg_V ] = pq_pv( pg,Z,Vp,V0 );
plot(pg(imag(qg_V)==0),qg_V(imag(qg_V)==0),'r','Linewidth',2);
axis([min(pg) max(pg) min(qg) max(qg)]);
xs = axis; axis equal;

plot(-real(S0),-imag(S0),'kx')
plot(real(Snb)*[1 1],xs(3:4),'k:')
plot(real(Snd)*[1 1],xs(3:4),'k-.')
plot(Pgh*[1 1],xs(3:4),'k--')
plot(real(Snp)*[1 1],xs(3:4),'k')

pg_curve = pg.*((pg<real(Snd)).*(pg>real(Snb)));
[ ~,qg_curve ] = pq_pv( pg_curve,Z,Vp,V0 );

Q_eps=1e-2;

plot([-real(S0),real(Snb)],...
     [-imag(S0)-2*Q_eps,imag(Snb)-2*Q_eps],'g','Linewidth',2);
plot([-real(S0),real(Snb),pg_curve(pg_curve>0)],...
     [-imag(S0)+Q_eps,imag(Snb)+Q_eps,qg_curve(pg_curve>0)],'b','Linewidth',2);
 
[ ~,~,~,Qm ] = sln_boundary( pg,Z,V0 );
plot(pg,Qm);

xlabel('$P_{n}$ (pu)'); ylabel('$Q_{n}$ (pu)');
lgnd = legend('$P_{t}$','$V_{+}$','$S_{0}$','$P_{n}^{nom}$',...
            '$\tilde{P}_{n}(\tilde{Q}_{g})$','$\hat{P}_{n}$','$P^{\prime}_{n}$',...
            '$S_{n}(\tau,0)$','$S_{n}(\tau, \tilde{Q}_{g})$');
set(lgnd,'FontSize',16,'Interpreter','Latex');

% export_fig(fig,fig_name);
% export_fig(fig,[fig_name,'.pdf'],'-dpdf');
