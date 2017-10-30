function [ idx,Pbar,Sdot,Shat,Sprm ] = find_qsln_idx_v2( PP,QQ,LL,VG,Qp,Vp,Ps,plt )

VpNaN_outs = 0./(VG<Vp) + 0./(VG>0.9); %avoid high and low voltages

if Qp==0
    QNaN_outs = zeros(size(QQ));
    QNaN_outs(1:end-1,:) = NaN;
    
    PP0 = PP + QNaN_outs + VpNaN_outs;
    QQ0 = QQ + QNaN_outs + VpNaN_outs;
    
    PPT = PP0 - real(LL);
    PPt = PPT(end,:);
else
    QNaN_outs = 0./(QQ>Qp);
    
    PP0 = PP + QNaN_outs + VpNaN_outs;
    QQ0 = QQ + QNaN_outs + VpNaN_outs;
    
    PPT = PP0 - real(LL);
    PPt = max(PPT,[],1);
end



[ Pbar,Pprm ] = find_pb_pp( PP,LL,VG,Vp );
Qprm = QQ0( PP0==Pprm );

PP0_ln = PP0(PPT==PPt);
PP0_ln = PP0_ln + 0./(PP0_ln<Pprm); 

if isinf(Qp)
    Pabs = abs(PP0_ln - max(Ps));
    Phat = PP0_ln(Pabs==min(min(Pabs)));
    Qhat = QQ0( PP0==Phat );
else
    Phat = max(Ps);
    Qhat = NaN;
end


Pdot = min( [max(PP0_ln);Phat] );
Qdot = QQ0( PP0==Phat );

Sdot = Pdot + 1i*Qdot;
Shat = Phat + 1i*Qhat;
Sprm = Phat + 1i*Qprm;

idx_0 = find(PP==PP(end,1));
idx_prm = find(PP0==Pprm);

idx = zeros(size(Ps));
for i = 1:numel(Ps)
    if Ps(i) == 0
        idx(i) = idx_0;
    elseif Ps(i) > Pprm
        idx(i) = idx_prm;
    else 
        Pabs = abs(PP0_ln - Ps(i));
        Pp_ln = PP0_ln(Pabs==min(min(Pabs)));
        if numel(Pp_ln)~=1
            QQQ = 1;
        end
        idxs = find( PP0==Pp_ln );
        if numel(idxs)==1
            idx(i) = idxs;
        else
            idx(i) = find( PPT==max(PPT(idxs)) );
        end
    end
end


if plt
    figure('Color','White','Position',[100 20 800 1000]);
    subplot(3,2,1:2);
    contour(PP,QQ,VG,'k'); axis equal; hold on;
    cc = contourf(PP0,QQ0,VG); clabel(cc);
    xs = axis;
    
    plot(Pbar*[1 1],xs(3:4),'k:');
    plot(Pdot*[1 1],xs(3:4),'k-.');
    plot(Phat*[1 1],xs(3:4),'k--');
    plot(Pprm*[1 1],xs(3:4),'k');
    
    subplot(323);
    plot(PP0_ln); hold on;
    plot(PPt);
    
    subplot(324);
    plot(Ps); hold on; plot(PP(idx)); hold on; plot(PPT(idx));
    subplot(325);
    plot(QQ(idx));
end

% [Pbar Pdot Phat Pprm]


end

