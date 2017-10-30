function [ Pb,Pp ] = find_pb_pp( PP,LL,VG,Vp )

VNaN_outs = 0./(VG<Vp);
PP = PP+VNaN_outs;
PT = PP-real(LL);

Pt_mx = max(max(PT));

Pp = PP(Pt_mx==PT);

Pb = max(PP(end,:));

end

