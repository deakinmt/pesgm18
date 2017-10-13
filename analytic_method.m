clear all; close all; clc
Z = exp(1i*pi/4);
V0 = 1.0;
S0 = 0.08+0.02*1i;
Vp = 1.06;
Qgen = 0.03;



pg = (-0.1:0.01:1.2);
qg = (-0.35:0.01:0.05);

[Pg,Qg] = meshgrid(pg,qg);

Sg = Pg + 1i*Qg;

[V2,Sl,~,~] = V2_Sl_calc( Sg,Z,V0,'p' );

[ Se,Sg,Snb,Snd,Snp ] = calc_s1( Sg,S0,Qgen(1),Z,Vp,V0 );

% contour(Pg,Qg,sqrt(V2),(0.6:0.05:1.1)); axis equal; hold on;
contour(Pg,Qg,real(Sl)); axis equal; hold on;
contour(Pg,Qg,sqrt(V2),(Vp:0.001:Vp+0.001),'r','Linewidth',2);
contour(Pg,Qg,sqrt(V2),(0.9:0.001:0.9+0.001),'b','Linewidth',2);
plot(-real(S0),-imag(S0),'kx')
plot(real(Snb),imag(Snb),'cx')
plot(real(Snd),imag(Snd),'gx')
plot(real(Snp),imag(Snp),'bx')

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







