function [ BXp,BXq ] = BXpq( XX,Z,Vg,V0 )

R = real(Z);
X = imag(Z);
aZ = abs(Z);


% BX = (X.*(Vg.^2)).^2 - (aZ.^2).*( ...
%                     (aZ.*Pnb).^2 - 2*Pnb.*R.*(Vg.^2) + (Vg.^2).*(Vg.^2 - V0.^2)...
%                                                     );
BXp = (R.*(Vg.^2)).^2 - (aZ.^2).*( ...
                    (aZ.*XX).^2 - 2*XX.*X.*(Vg.^2) + (Vg.^2).*(Vg.^2 - V0.^2)...
                                                    );

BXq = (X.*(Vg.^2)).^2 - (aZ.^2).*( ...
                    (aZ.*XX).^2 - 2*XX.*R.*(Vg.^2) + (Vg.^2).*(Vg.^2 - V0.^2)...
                                                    );

end

