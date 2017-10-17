function [ Pp,Pm,Qp,Qm ] = sln_boundary( xx,Z,V0 )

R = real(Z);
X = imag(Z);

aP = X^2;
bP = -(R*(V0^2) + 2*xx*R*X);
cP = (xx.^2)*(R^2) - ( (V0^4)/4 ) - (V0^2)*xx*X ;

aQ = R^2;
bQ = -(X*(V0^2) + 2*xx*R*X);
cQ =  (xx.^2)*(X^2) - ( (V0^4)/4 ) - (V0^2)*xx*R;

DP = ( (bP.^2) - 4*aP.*cP );
DQ = ( (bQ.^2) - 4*aQ.*cQ );

Pp = (-bP + sqrt(DP))./(2*aP);
Pm = (-bP - sqrt(DP))./(2*aP);
Qp = (-bQ + sqrt(DQ))./(2*aQ);
Qm = (-bQ - sqrt(DQ))./(2*aQ);

end

