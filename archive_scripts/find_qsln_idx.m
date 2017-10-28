function [ idx,Sbar,Sdot,Shat ] = find_qsln_idx( PP,QQ,VV,VG,Qp,Vp,Ps,dQ,dV,plt )

% first find the unconstrained reactive power:
VpNaN_outs = 0./(VG<Vp);
VloNaN_outs = 0./(VG>(Vp-dV));
% VpNaN_outs = 0./(VV<Vp);
% VloNaN_outs = 0./(VV>(Vp-dV));
PNaN_outs = 0./(abs(PP)<max(Ps));


PVhatmat  = PP + PNaN_outs + VpNaN_outs + VloNaN_outs;
Phat = PP(PP==max(max(PVhatmat)));
Qhat = QQ(PP==max(max(PVhatmat)));

% now find the output real/reactive powers for a given Ps
Qp = min(abs([Qp,Qhat]));
NoQNaN_outs = 0./(abs(QQ)<dQ);
QNaN_outs = 0./(abs(QQ)<Qp);
% spy((VmaxMat>(Vp-dV)).*(VmaxMat<Vp))

Pgnsetmat = PP + QNaN_outs + VpNaN_outs;
Qgnsetmat = QQ + QNaN_outs + VpNaN_outs;
Sgnsetmat = Pgnsetmat + 1i*Qgnsetmat;

PVlinemat = PP + QNaN_outs + VpNaN_outs + VloNaN_outs;
QVlinemat = QQ + QNaN_outs + VpNaN_outs + VloNaN_outs;
PNoQgnmat = PP + QNaN_outs + VpNaN_outs + NoQNaN_outs;
QNoQgnmat = QQ + QNaN_outs + VpNaN_outs + NoQNaN_outs;

Pdot = PP(PP==max(max(Pgnsetmat)));
Qdot = QQ(PP==max(max(Pgnsetmat)));
Pbar = PP(PP==max(max(PNoQgnmat)));
Qbar = QQ(PP==max(max(PNoQgnmat)));

idx = zeros(size(Ps));
for i = 1:numel(Ps)
    if Ps(i) == 0
        Sabs = abs(Sgnsetmat);
        idx(i) = find( Sabs==min(min(Sabs)) );
    elseif Ps(i) < Pbar
        Pabsb = abs(PNoQgnmat-Ps(i));
        idx(i) = find( Pabsb==min(min(Pabsb)) );
    elseif Ps(i) < Pdot
        Pabsd = abs(PVlinemat-Ps(i));
        idx(i) = find( Pabsd== min(min(Pabsd)) );
    else
        idx(i) = find( PP==Pdot(1),1 );
    end
end


if plt
    figure;
    contour(PP,QQ,VG,'k'); axis equal; hold on;
    cc = contour(Pgnsetmat,Qgnsetmat,VG); clabel(cc);
    cc = contourf(PVlinemat,QVlinemat,VG); clabel(cc);
    cc = contourf(PNoQgnmat,QNoQgnmat,VG); clabel(cc);

    xs = axis;
    plot(Pbar*[1 1],xs(3:4),'k:');
    plot(Pdot*[1 1],xs(3:4),'k-.');
    plot(Phat*[1 1],xs(3:4),'k--');
end

Sbar = Pbar + 1i*Qbar;
Sdot = Pdot + 1i*Qdot;
Shat = Phat + 1i*Qhat;
end

