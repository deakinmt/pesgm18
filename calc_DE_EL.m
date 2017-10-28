function [ Qn,D_En,D_Et,e_L ] = calc_DE_EL( Ps0,PPg,QQg,Vmx,Vgn,Vp,n,LL )

% inputs
% Ps0 irradiance time series
% PP real power generation matrix
% QQ reactive power generation matrix
% Vmx maximum voltage on feeder
% Vgn generator voltage matrix
% Vp voltage limit
% LL (complex) loss matrix

dQ = 2e-2;
dV = 8e-3;
dt = 1/60; % assume one minute resolution

[ idx{1},~,~,~ ] = find_qsln_idx( PPg,QQg,Vmx,Vgn,dQ,Vp,Ps0,dQ,dV,0 );
[ idx{n+1},~,~,ShatI ] = find_qsln_idx( PPg,QQg,Vmx,Vgn,inf,Vp,Ps0,dQ,dV,1 );
QhatI = imag(ShatI);

PPt = PPg - real(LL);

Qn = zeros(n,1);
Qn(end) = QhatI;

for i = 2:n
    Qn(i-1) = QhatI*((i-1)/n);
    [ idx{i},~,~,~ ] = find_qsln_idx( PPg,QQg,Vmx,Vgn,Qn(i-1),Vp,Ps0,dQ,dV,0 );
end

D_En = zeros(n,1);
D_Et = zeros(n,1);

for i = 1:n
    D_En(i) = dt*sum( PPg(idx{i}) - PPg(idx{1}) );
    D_Et(i) = dt*sum( PPt(idx{i}) - PPt(idx{1}) );
end

e_L = 100*(D_En-D_Et)./D_Et;




end

