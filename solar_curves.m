% solar curves
close all; clear all; clc
fig_loc = 'C:\Users\chri3793\Documents\DPhil\pesgm18\pesgm18_paper\figures\';

set(0,'defaulttextinterpreter','latex');
set(0,'defaultaxesfontsize',14);
fig_nompos = [100 100 550 320];

%% First create some solar data:

% delta = 23.65*pi/180; % july day
% 
% atan(cot(delta)/1e-3)*180/pi % North pole
% atan(cot(delta)/1)*180/pi % Lulea, Fairbanks, Reykjavik
% atan(cot(delta)/2)*180/pi % Paris, Winnipeg
% atan(cot(delta)/4)*180/pi % Houston, New Dehli
% atan(cot(delta)/1e3)*180/pi % Nairobi, Singapore
t = (-12:1/(60):12-1/(60));
w = pi/12; %rad/hour

k = 2;
Ps_pk = 1.0;

Ps = Ps_pk*(1 + k*cos(w*t))/(1+k);
Ps = Ps.*(Ps>0);
% save('solar_curves_dataset.mat','Ps','t');
%%
Pgp = 1.08;
Pgh = 1.0;
Pgd = 0.78;
Pgb = 0.55;

Ps_b = [Ps; Pgb*ones(size(Ps))];
Pshade00 = zeros(size(Ps))./(Ps~=0);
Pshade01 = min( [Ps;Pgb*ones(size(Ps))] ,[],1)./(Ps~=0);
Pshade10 = (Pgb)*ones(size(Ps))./(Ps>Pgb);
Pshade11 = min( [Ps;Pgd*ones(size(Ps))] ,[],1)./(Ps>Pgb);
% Pshade10 = Pgd*ones(size(Ps))./(Ps>Pgd);
% Pshade11 = Ps./(Ps>Pgd);

fig_name = [fig_loc,'solar_curves'];
fig = figure('Color','White','Position',fig_nompos);

% p4 = plot(t+12,Ps,'Linewidth',2,'Color',[0.6 0.3 0.5]); hold on;
p6=shadedplot(t+12,Pshade10,Pshade11,[0.8 0.8 0.98],'w'); hold on;
p5=shadedplot(t+12,Pshade00,Pshade01,[0.8 0.98 0.8],'g'); hold on;
xs = axis;

p0 = plot(xs(1:2),Pgp*[1 1],'k');
p1 = plot(xs(1:2),Pgh*[1 1],'k--');
p2 = plot(xs(1:2),Pgd*[1 1],'k-.');
p3 = plot(xs(1:2),Pgb*[1 1],'k:');

% arrow([3 Pgd],[3 Pgd-0.1]);
% arrow([3 Pgd],[3 Pgd+0.1]);
p4 = plot(t+12,Ps,'Linewidth',2,'Color',[0.6 0.3 0.5]); hold on;

xlabel('$\tau $ (hour)');
ylabel('$P_{g} / \hat{P}_{g}$');
lgnd = legend([p4,p5(2),p6(2),p0,p1,p2,p3],'$P_{S}(\tau)$','$\bar{E}_{g}$','$\Delta E_{g}$','$P^{\prime}_{g}$',...
            '$\hat{P}_{g}$','$\dot{P}_{g}(\dot{Q}_{g})$','$\bar{P}_{g}$');
set(lgnd,'Interpreter','Latex','FontSize',14);
axis([0 24 0 1.3]);

xticks([0 6 12 18 24]);
grid off;
% export_fig(fig,fig_name);
% export_fig(fig,[fig_name,'.pdf'],'-dpdf');








