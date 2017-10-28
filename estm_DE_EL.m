function [ Qn,D_En,D_Et,e_L ] = estm_DE_EL( Ps0,Ps_k,S0,Z,Vp,V0,n )


PsMax = (max(Ps0)*Ps_k) - real(S0);

[ ~,Qnd ] = pq_pv( PsMax,Z,Vp,V0 );
Qn = linspace(1/n,1,n)'*Qnd;

dt = 1/60; %assume one minute data

% St = zeros(numel(Ps0),numel(Ps_k));
% Sg = zeros(numel(Ps0),numel(Ps_k));

D_En = zeros(n,numel(Ps_k));
D_Et = zeros(n,numel(Ps_k));

for i = 1:n
    for j = 1:numel(Ps_k)
        
        if i==1
            [ St0,Sg0,~,~,~ ] = calc_s1( Ps0*Ps_k(j),S0,0,Z,Vp,V0 );
            Pg0 = sum(real(Sg0));
            Pt0 = sum(real(St0));
        end
        
        [ St,Sg,~,~,~ ] = calc_s1( Ps0*Ps_k(j),S0,-Qn(i,j),Z,Vp,V0 );
        D_En(i,j) = dt*(sum(real(Sg)) - Pg0);
        D_Et(i,j) = dt*(sum(real(St)) - Pt0);
    end
end

e_L = 100*(D_En - D_Et)./D_Et;


