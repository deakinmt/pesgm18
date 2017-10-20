clc

V0 = [0.9 1.05 1.1];
Vp = [1.1 1.06  0.9];
lv = Vp./V0;


lz_p0_0 = (V0.^2)./sqrt( 4*(Vp.^2) - (V0.^2) )
lz_p0= V0./sqrt( 4*(Vp.^2) - (V0.^2) )
% lz_p = 1./sqrt( 4*(lv.^2) - 1)

% lz_c = lv./sqrt( 4 - (lv.^2) )
% 
% lr_c = lz_c./sqrt(1 + lz_c.^2)  
% lr_p = lz_p./sqrt(1 + lz_p.^2)
% 
% lr_c.*lr_p

%%
lz_i = sqrt( lv.^2 - 1 )



