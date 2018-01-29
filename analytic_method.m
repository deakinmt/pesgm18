clear all; close all; clc

addpath('pesgm_funcs');

% fig_loc = 'C:\Users\Matt\Documents\DPhil\pesgm18\pesgm18_paper\figures\';
fig_loc = 'C:\Users\chri3793\Documents\DPhil\pesgm18\pesgm18_paper\figures\';

set(0,'defaulttextinterpreter','latex');
set(0,'defaultaxesfontsize',14);
set(0,'defaulttextfontsize',14);

fig_nompos = [100 100 550 300];

%% 
Vp = 1.06;
V0 = 1.0;
Z = 0.203*exp(1i*acot(1.85)); %see PSCC paper
S0 = 0.72*exp(1i*9.4*pi/180);
% S0 = 0;
Qgen = 0.8;

ds = 0.003;

pg = 5*(-0.2:ds:0.75);
% qg = 5*(-0.6:ds:0.15);
qg = 5*(-0.6:ds:0.15);

[Pg,Qg] = meshgrid(pg,qg);

Sg = Pg + 1i*Qg;

[~,Sl,~,~] = V2_Sl_calc( Sg,Z,V0,'p' );
[ Se,Sg,Snb,Snd,Snp ] = calc_s1( Sg,S0,Qgen(1),Z,Vp,V0 );
Pgh = 0.8*real(Snp);
[~,Qgh] = pq_pv( Pgh,Z,Vp,V0 );
%
%------------------ PLOT
fig_name = [fig_loc,'analytic_method'];
fig = figure('Color','White');
[~,PP] = contour(Pg,Qg,Pg-real(Sl)); axis equal; hold on;
[ ~,qg_V ] = pq_pv( pg,Z,Vp,V0 );
% PP(2) = plot(pg(imag(qg_V)==0),qg_V(imag(qg_V)==0),'r','Linewidth',2);
PP(2) = plot(pg(imag(qg_V)==0),qg_V(imag(qg_V)==0),'r');
axis([min(pg) max(pg) min(qg) max(qg)]);
xs = axis; axis equal;

pg_curve = pg.*((pg<real(Snd)).*(pg>real(Snb)));
[ ~,qg_curve ] = pq_pv( pg_curve,Z,Vp,V0 );

Q_eps=1e-2;

PP(3) = plot([-real(S0),real(Snb)],...
     [-imag(S0)-2*Q_eps,imag(Snb)-2*Q_eps],'g--','Linewidth',2);
PP(4) = plot([-real(S0),real(Snb),pg_curve(pg_curve>0)],...
     [-imag(S0)+Q_eps,imag(Snb)+Q_eps,qg_curve(pg_curve>0)],'b-.','Linewidth',2);

ds=0.04;

plot(-real(S0),-imag(S0),'ko'); text(-real(S0)+ds,-imag(S0)+3.5*ds,'$S_{0}$');
plot(real(Snb),imag(Snb),'k+'); text(real(Snb)+ds,imag(Snb)+2*ds,'${S}_{n}$\textsuperscript{nom}');
plot(real(Snd),imag(Snd),'k+'); text(real(Snd)+ds,imag(Snd)+2*ds,'$\tilde{S}_{n}(\tilde{Q}_{n})$');
plot(Pgh,Qgh,'k+'); text(Pgh+ds,Qgh+ds+2*ds,'$\hat{S}_{n}$');
plot(real(Snp),imag(Snp),'k+'); text(real(Snp)+ds,imag(Snp)+2*ds,'$S^{\prime}_{n}$');

% PP(5) = plot(real(Snb)*[1 1],xs(3:4),'k:');
% PP(6) = plot(real(Snd)*[1 1],xs(3:4),'k-.');
plot(real(Snb)*[1 1],xs(3:4),'k:');
plot(real(Snd)*[1 1],xs(3:4),'k-.');
% plot(Pgh*[1 1],xs(3:4),'k--')
% plot(real(Snp)*[1 1],xs(3:4),'k')


[ ~,~,~,Qm ] = sln_boundary( pg,Z,V0 );
plot(pg,Qm,'k');

xlabel('Net Real Power, $P_{n}$ (pu)'); ylabel('Net Reactive Power, $Q_{n}$ (pu)');
% lgnd = legend(PP,'$P_{t}$','$V_{+}$','$S_{n}(\tau,0)$','$S_{n}(\tau, \tilde{Q}_{g})$','$\bar{P}_{n}$','$\tilde{P}_{n}(\tilde{Q}_{n})$');
lgnd = legend(PP,'$P_{t}$','$V_{+}$','$S_{n}(\tau,0)$','$S_{n}(\tau, \tilde{Q}_{g})$');
set(lgnd,'FontSize',16,'Interpreter','Latex');

text(0.3,-2.95,'${P}_{n}$\textsuperscript{nom}','Interpreter','Latex','Rotation',90);
text(0.86,-2.95,'$\tilde{P}_{n}(\tilde{Q}_{n})$','Interpreter','Latex','Rotation',90);

% export_fig(fig,fig_name);
% export_fig(fig,[fig_name,'.pdf'],'-dpdf');

%% OLD VERSION
Vp = 1.06; V0 = 1.0;
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

%------------------ PLOT ===> %% OLD VERSION
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


