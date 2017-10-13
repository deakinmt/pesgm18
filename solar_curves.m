% solar curves
clc
delta = 23.65*pi/180;

atan(cot(delta)/1)*180/pi % Lulea, Fairbanks, Reykjavik
atan(cot(delta)/2)*180/pi % Paris, Winnipeg
atan(cot(delta)/4)*180/pi % Houston, New Dehli
atan(cot(delta)/1e3)*180/pi % Nairobi, Singapore

%%
X0 = (0:1:90);
k = 2;
plot(sin(X0*pi/180).*sin(23.65*pi/180)); hold on;
plot(cos(X0*pi/180).*cos(23.65*pi/180));

%%
x0 = 49; %e.g. paris, bavaria, and a large array (e.g. on a rooftop).


sin(x0*pi/180).*sin(23.65*pi/180);
cos(x0*pi/180).*cos(23.65*pi/180);





x = (-2*pi/3:0.01:2*pi/3);
plot(0.33 + 0.66*cos(x)); % 40% capacity factor






