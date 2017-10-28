function [ St,Sg,Snb,Snd,Snp ] = calc_s1( Ps,S0,Qgen,Z,Vp,V0 )

P0 = real(S0); Q0 = imag(S0);


[ Snb,Snd,Snp ] = calc_xpts( Q0,Qgen,Z,Vp,V0 );

Pnb = real(Snb); Pnd = real(Snd); Qnd = imag(Snd);

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
St = Sn - Sl;


end

