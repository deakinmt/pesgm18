function [ Qgd,D_Eg,D_Et,e_L ] = estm_DE_EL( Ps0,Ps_k,S0,Z,Vp,V0,n )

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

% for i = 1:n
%     for j = 1:numel(Ps_k)
%         Ps = Ps0*( Pgb + Ps_k(j)*( Pgp - Pgb ) );
%         if i==1
%             [ St0,Sg0,~,~,~ ] = calc_s1( Ps,S0,0,Z,Vp,V0 );
%             Pg0 = sum(real(Sg0));
%             Pt0 = sum(real(St0));
%         end
%         
%         [ St,Sg,~,~,~ ] = calc_s1( Ps,S0,-Qgd(i,j),Z,Vp,V0 );
%         D_Eg(i,j) = dt*(sum(real(Sg)) - Pg0);
%         D_Et(i,j) = dt*(sum(real(St)) - Pt0);
%     end
% end
% 
% e_L = 100*(D_Eg - D_Et)./D_Et;


