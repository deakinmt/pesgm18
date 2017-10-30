function [ Qgd,D_Eg,D_Et,e_L ] = calc_DE_EL( Ps0,PPg,QQg,Vmx,Vgn,Vp,n,LL )

% inputs
% Ps0 irradiance time series
% PP real power generation matrix
% QQ reactive power generation matrix
% Vmx maximum voltage on feeder
% Vgn generator voltage matrix
% Vp voltage limit
% LL (complex) loss matrix
dt = 1/60; % assume one minute resolution, going to MWh.

[ idx{1},~,~,~ ] = find_qsln_idx_v2( PPg,QQg,LL,Vgn,0,Vp,Ps0,0 );
[ idx{n+1},~,~,ShatI ] = find_qsln_idx_v2( PPg,QQg,LL,Vgn,-inf,Vp,Ps0,0 );

QhatI = imag(ShatI);

Qgd = linspace(1/n,1,n)'*QhatI;
for i = 1:n-1
    [ idx{i+1},~,~,~ ] = find_qsln_idx_v2( PPg,QQg,LL,Vgn,Qgd(i),Vp,Ps0,0 );
end

D_Eg = zeros(n,1);
D_Et = zeros(n,1);

PPt = PPg - real(LL); %pu

PPg0 = dt*sum( PPg(idx{1}) ); %pu.hour
PPt0 = dt*sum( PPt(idx{1}) ); %pu.hour
for i = 1:n
    D_Eg(i) = dt*sum( PPg(idx{i+1}) ) - PPg0; %pu.hour
    D_Et(i) = dt*sum( PPt(idx{i+1}) ) - PPt0; %pu.hour
end

e_L = 100*(D_Eg-D_Et)./D_Et;




end

