function [ Se,Sg,Snb,Snd,Snp ] = calc_s1( Ps,S0,Qgen,Z,Vp,V0 )

P0 = real(S0); Q0 = imag(S0);

% first calculate the crossing points
[ Pnb,~ ] = pq_pv( -Q0,Z,Vp,V0 );
Snb = Pnb - 1i*Q0;
[ Pnd,~ ] = pq_pv( - Q0 - Qgen,Z,Vp,V0 );
Qnd = -(Q0 + Qgen);
Snd = Pnd - 1i*(Q0 + Qgen);
[ Pnp,Qnp,~ ] = pscc18_theorem_1( Vp, V0, Z );
Snp = Pnp + 1i*Qnp;


% now calculate the exported powers. This involves calculating the
% curtailed power.
Pg = zeros(size(Ps)); Qg = zeros(size(Ps));

for i = 1:numel(Ps)
    if Ps(i)-P0 <= Pnb
        Pg(i) = Ps(i);
        Qg(i) = 0;
        
    elseif Ps(i)-P0 <= Pnd
        Pg(i) = Ps(i);
        [ ~,Qn ] = pq_pv( Pg(i) - P0,Z,Vp,V0 );
        Qg(i) = Qn+Q0;
        
    elseif Ps(i)-P0 > Pnd
        Pg(i) = Pnd + P0;
        Qg(i) = Qnd + Q0;
    end
end

Sg = Pg + 1i*Qg;
Sn = Sg - S0;
Pc = Ps - Pg;
[~,Sl,~,~] = V2_Sl_calc( Sn,Z,V0,'p' );
Se = Sn - Sl;


end

