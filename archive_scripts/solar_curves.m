% solar curves
clc
fig_loc = 'C:\Users\chri3793\Documents\DPhil\malcolm_updates\wc171016\figures\';

delta = 23.65*pi/180;

atan(cot(delta)/1e-3)*180/pi % North pole
atan(cot(delta)/1)*180/pi % Lulea, Fairbanks, Reykjavik
atan(cot(delta)/2)*180/pi % Paris, Winnipeg
atan(cot(delta)/4)*180/pi % Houston, New Dehli
atan(cot(delta)/1e3)*180/pi % Nairobi, Singapore

%%
fig_name = [fig_loc,'solar_curves'];
fig = figure('Color','White');
t = (-12:1/(60):12-1/(60));
w = pi/12; %rad/hour

k = 2;
Ps_pk = 1.0;

Pgd = 0.55;

Ps = Ps_pk*(1 + k*cos(w*t))/(1+k);
Ps = Ps.*(Ps>0);
Pshade1 = Pgd*ones(size(Ps))./(Ps>Pgd);
Pshade2 = Ps./(Ps>Pgd);

plot(t+12,Ps,'Linewidth',2); hold on;

t1 = acos( (Pgd*(1 + k) - 1)/k )/w;
xs = axis;

plot(xs(1:2),Pgd*[1 1],'k--');
plot(t1*[1 1]+12,xs(3:4),'k:');
plot(-t1*[1 1]+12,xs(3:4),'k:');

arrow([3 Pgd],[3 Pgd-0.1]);
arrow([3 Pgd],[3 Pgd+0.1]);
shadedplot(t+12,Pshade1,Pshade2,[0.5 0.3 0.5],'k'); hold on;

xlabel('Hour');
ylabel('Pn/Hat(Pn)');
legend('P_S','Dot(Pn)','t1'); grid on;
save('solar_curves_dataset.mat','Ps','t');
% export_fig(fig,fig_name);
% export_fig(fig,[fig_name,'.pdf'],'-dpdf');
%%
fig_name = [fig_loc,'solar_curves_curt'];

Pgb = 0.6;

PGD = (0:0.001:1.1);
T1 = acos( ( PGD*(1+k) - 1 )/k )/w;

% Ecurt = (2*k/(1+k))*( (sin(w*T1)/w) - T1.*cos(w*T1) );

% plot(PGD,gradient(Ecurt));
plot(PGD,real(T1));

dEdPGD = -(2/w)*T1;
% plot(PGD,dEdPGD);

%%
x0 = 49; %e.g. paris, bavaria, and a large array (e.g. on a rooftop).
t = (-12:0.05:11.5); 

sin(x0*pi/180).*sin(23.65*pi/180);
cos(x0*pi/180).*cos(23.65*pi/180);
t = (-2*pi/3:0.01:2*pi/3);
plot(0.33 + 0.66*cos(t)); % 40% capacity factor






