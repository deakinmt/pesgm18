function [ Snb,Snd,Snp ] = calc_xpts( Q0,Qg,Z,Vp,V0 )

[ Pnb,~ ] = pq_pv( -Q0,Z,Vp,V0 );
Snb = Pnb - 1i*Q0;
[ Pnd,~ ] = pq_pv( - Q0 - Qg,Z,Vp,V0 );
Qnd = -(Q0 + Qg);
Snd = Pnd - 1i*(Q0 + Qg);
[ Pnp,Qnp,~ ] = pscc18_theorem_1( Vp, V0, Z );
Snp = Pnp + 1i*Qnp;




end

