function [ ppv,qpv ] = pq_pv( xx,Z,Vg,V0 )

[ ~,lr,lx,~,~ ] = lambdas( Z );
lv = Vg./V0;
Sbar = (V0.^2)./abs(Z);


ppv = Sbar.*( (lr.*(lv.^2)) - sqrt( lv.^2 - (lv.^4.*lx.^2) + (2*lv.^2.*lx.*xx./Sbar) ...
                                                        - (xx./Sbar).^2 ) );

qpv = Sbar.*( (lx.*(lv.^2)) - sqrt( lv.^2 - (lv.^4.*lr.^2) + (2*lv.^2.*lr.*xx./Sbar) ...
                                                        - (xx./Sbar).^2 ) );



end