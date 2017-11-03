function [ Qgd,D_Eg,D_Et,e_L,eps_P,DPg ] = estm_DE_EL( Ps0,Ps_k,S0,Z,Vp,V0,n )

[ Pnp,Qnp,~ ] = pscc18_theorem_1( Vp, V0, Z );
[ Snb,~,~ ] = calc_xpts( imag(S0),Qnp,Z,Vp,V0 );
Pnb = real(Snb);

Pgp = Pnp + real(S0);
Pgb = Pnb + real(S0);
% Qgp = Qnp + imag(S0);

PnMax = (Pgb + Ps_k*(Pgp-Pgb)) - real(S0) ;
[ ~,Qnd ] = pq_pv( PnMax ,Z,Vp,V0 );

Qgd = linspace(1/n,1,n)'*( Qnd + imag(S0) );

dt = 1/60; %assume one minute data

% St = zeros(numel(Ps0),numel(Ps_k));
% Sg = zeros(numel(Ps0),numel(Ps_k));

D_Eg = zeros(n,numel(Ps_k));
D_Et = zeros(n,numel(Ps_k));

for i=1:numel(Ps_k)
    Ps = Ps0*( Pgb + Ps_k(i)*( Pgp - Pgb ) );
    
    [ St0,Sg0,~,~,~ ] = calc_s1( Ps,S0,0,Z,Vp,V0 );
    Pg0 = sum(real(Sg0));
    Pt0 = sum(real(St0));

    for j = 1:n
        [ St,Sg,~,~,~ ] = calc_s1( Ps,S0,-Qgd(j,i),Z,Vp,V0 );
        D_Eg(j,i) = dt*(sum(real(Sg)) - Pg0);
        D_Et(j,i) = dt*(sum(real(St)) - Pt0);
    end
    
end

e_L = 100*(D_Eg - D_Et)./D_Et;



%% calculate an estimate of the utility power error

QNd = Qgd - imag(S0);
% QNd = linspace(1/n,1,n)'*( Qnd );
[ PNd,~ ] = pq_pv( QNd ,Z,Vp,V0 );
[~,~,Pl,~] = V2_Sl_calc( PNd+1i*QNd,Z,V0,'p' );

DPg = PNd - PNd(1);
DPl = Pl - Pl(1);
DPt = DPg - DPl;
eps_P = (DPg - DPt)./DPt;

% plot(Qgd/min(Qgd),eps_P);



