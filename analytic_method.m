clear all; close all; clc

fig_loc = 'C:\Users\chri3793\Documents\DPhil\malcolm_updates\wc171016\figures\';

Vp = 1.00;
V0 = 0.9;
Z = 0.203*exp(1i*acot(1.85)); %see PSCC paper
% S0 = 0.72*exp(1i*9.4*pi/180);
S0 = 0;
Qgen = 1;

ds = 0.003;

pg = 5*(-0.2:ds:0.75);
qg = 5*(-0.6:ds:0.15);

[Pg,Qg] = meshgrid(pg,qg);

Sg = Pg + 1i*Qg;

[~,Sl,~,~] = V2_Sl_calc( Sg,Z,V0,'p' );
[ Se,Sg,Snb,Snd,Snp ] = calc_s1( Sg,S0,Qgen(1),Z,Vp,V0 );
Pgh = 0.8*real(Snp);

fig_name = [fig_loc,'analytic_method'];
fig = figure('Color','White');
contour(Pg,Qg,Pg-real(Sl)); axis equal; hold on;
% contour(Pg,Qg,Pg-real(Sl),(0.99:0.02:1.01),'b','Linewidth',2); 
[ ~,qg_V ] = pq_pv( pg,Z,Vp,V0 );
plot(pg(imag(qg_V)==0),qg_V(imag(qg_V)==0),'r','Linewidth',2);
axis([min(pg) max(pg) min(qg) max(qg)]);
xs = axis;

plot(-real(S0),-imag(S0),'kx')
plot(real(Snb)*[1 1],xs(3:4),'k:')
plot(real(Snd)*[1 1],xs(3:4),'k-.')
plot(Pgh*[1 1],xs(3:4),'k--')
plot(real(Snp)*[1 1],xs(3:4),'k')

% plot the direction curve
pg_curve = pg.*((pg<real(Snd)).*(pg>real(Snb)));
[ ~,qg_curve ] = pq_pv( pg_curve,Z,Vp,V0 );

plot([-real(S0),real(Snb),pg_curve(pg_curve>0)],...
     [-imag(S0),imag(Snb),qg_curve(pg_curve>0)],'g','Linewidth',2);
 
[ ~,~,~,Qm ] = sln_boundary( pg,Z,V0 );
plot(pg,Qm);

xlabel('Pn'); ylabel('Qn');
legend('Pt','Vp','S0','Pnb','Pnd','Pnh','Pnp','Straj');

% export_fig(fig,fig_name);
% export_fig(fig,[fig_name,'.pdf'],'-dpdf');
%%
Z = 0.203*exp(1i*acot(0.5)); %see PSCC paper
Qgen = 0.3;

ds = 0.003;

pg = 5*(-0.2:ds:1.2);
qg = 5*(-0.6:ds:0.15);

[Pg,Qg] = meshgrid(pg,qg);

Sg = Pg + 1i*Qg;

[~,Sl,~,~] = V2_Sl_calc( Sg,Z,V0,'p' );
[ Se,Sg,Snb,Snd,Snp ] = calc_s1( Sg,S0,Qgen(1),Z,Vp,V0 );
Pgh = 0.8*real(Snp);

fig_name = [fig_loc,'analytic_method_v2'];
fig = figure('Color','White');
contour(Pg,Qg,Pg-real(Sl)); axis equal; hold on;
[ ~,qg_V ] = pq_pv( pg,Z,Vp,V0 );
plot(pg(imag(qg_V)==0),qg_V(imag(qg_V)==0),'r','Linewidth',2);
axis([min(pg) max(pg) min(qg) max(qg)]);
xs = axis;

plot(-real(S0),-imag(S0),'kx')
plot(real(Snb)*[1 1],xs(3:4),'k:')
plot(real(Snd)*[1 1],xs(3:4),'k-.')
plot(Pgh*[1 1],xs(3:4),'k--')
plot(real(Snp)*[1 1],xs(3:4),'k')

% plot the direction curve
pg_curve = pg.*((pg<real(Snd)).*(pg>real(Snb)));
[ ~,qg_curve ] = pq_pv( pg_curve,Z,Vp,V0 );

plot([-real(S0),real(Snb),pg_curve(pg_curve>0)],...
     [-imag(S0),imag(Snb),qg_curve(pg_curve>0)],'g','Linewidth',2);
 
[ ~,~,~,Qm ] = sln_boundary( pg,Z,V0 );
plot(pg,Qm);

xlabel('Pn'); ylabel('Qn');
legend('Pt','Vp','S0','Pnb','Pnd','Pnh','Pnp','Straj','Location','SouthWest');


% export_fig(fig,fig_name);
% export_fig(fig,[fig_name,'.pdf'],'-dpdf');





%%
x = (-12:0.05:11.5); 
k = [0.5,1,2,1e3];
Qgen = (0:0.001:0.2);

Ps_pk = 0.4;
% Ps = Ps_pk*(1 + k*cos(x*pi/12))/(1+k);
% Ps = Ps.*(Ps>0);
% plot((0:0.5:23.5),Ps);
%
Eqgen = zeros(size(Qgen));
for i = 1:numel(k)
    Ps = Ps_pk*(1 + k(i)*cos(x*pi/12))/(1+k(i));
    Ps = Ps.*(Ps>0);
    for j = 1:numel(Qgen)
        [ Se,Sg,~,~,~ ] = calc_s1( Ps,S0,Qgen(j),Z,Vp,V0 );
        Eqgen(j) = sum(real(Se));
    end
    subplot(211)
    plot(Qgen,Eqgen - Eqgen(1)); hold on;
    subplot(212)
    plot(Qgen(2:end),Eqgen(2:end) - Eqgen(1:end-1)); hold on;
end

%%
sum(real(Se));

Sn = Sg - S0;
Pc = real(Sg) - Ps;
plot(real(Se)); hold on; plot(real(Sn)); grid on;
plot(imag(Sg)); plot(Pc);

legend('Se','Sn','Qg','Pc');







