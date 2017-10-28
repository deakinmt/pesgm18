function [ RR,figs ] = run_pesgm_feeder( FF,pl_options )

Vp = FF.Vp;
fig_nompos = FF.fig_nompos;

%% Find nominal network parameters from opendss

[~, DSSObj, DSSText] = DSSStartup;
DSSCircuit=DSSObj.ActiveCircuit; DSSSolution=DSSCircuit.Solution;
DSSText.command=['Compile (',pwd,FF.filename,'.dss)'];
DSSText.command=['edit vsource.source pu=',num2str(FF.V0)];
DSSSolution.Solve;

Zbus = create_zbus(DSSCircuit);
V0 = DSSCircuit.Vsource.pu;

DSSCircuit.SetActiveBus(FF.NUT);
vbN = sqrt(3)*DSSCircuit.ActiveBus.kVBase; % LL (@ NUT)
DSSCircuit.SetActiveBus(FF.SRC);
vbS = sqrt(3)*DSSCircuit.ActiveBus.kVbase; % LL (@ SRC)
nTr0 = vbN/vbS;
DSSCircuit.Transformers.First;
sb = DSSCircuit.Transformers.kVA/10; %ieee34mod1 script is factor of 10 out
zbN = 1e3*(vbN^2)/sb; %ohms
ibN = sb*sqrt(1/3)/vbN; %A, per phase

Z_ohm  = find_node_Z_pscc( DSSCircuit.YNodeOrder,Zbus,FF.SRC,FF.NUT,nTr0 );
Z = Z_ohm/zbN; %pu
Ssc = ((vbN^2)/abs(Z*zbN))*1e3; %kW
Sload = (1769 + 294*1i)/sb; %feeder nominal load power

% create generator at given node:
DSSText.Command=['new generator.gen phases=3 bus1=',FF.NUT,'.1.2.3 kw=0 pf=0'];
DSSText.Command='new monitor.genpq element=generator.gen terminal=1 mode=65 PPolar=no';
DSSText.Command='new monitor.genvi element=generator.gen terminal=1 mode=32 VIPolar=yes';
%% Run the DSS
[Pg, Qg] = meshgrid(Ssc*FF.pg_ssc,Ssc*FF.qg_ssc);
Pf = sign(Qg).*Pg./sqrt(Pg.^2 + Qg.^2);

totlss = zeros(numel(Pg),2);
totpwr = zeros(numel(Pg),2);
VgenMes = zeros(numel(Pg),1); 
VmaxMes = zeros(numel(Pg),1); 
DSSCircuit.Monitors.ResetAll;

aa = exp(1i*pi*2/3);
AA = [1;aa^2;aa]/sqrt(3);
NBus = find_node_idx(DSSCircuit.YNodeOrder,FF.NUT);

tic
for i = 1:numel(Pg)
    DSSText.Command=strcat('edit generator.gen kw=',num2str(Pg(i)));
    DSSText.Command=strcat('edit generator.gen pf=',num2str(Pf(i)));
    
    DSSSolution.Solve;
    DSSCircuit.Sample;
    
    
    YNodeVarray = DSSCircuit.YNodeVarray';
    YNodeV = YNodeVarray(1:2:end) + 1i*YNodeVarray(2:2:end);
    VgenMes(i) = abs(AA'*YNodeV(NBus))/(vbN*1e3);
    
    AllBusVmagPu = DSSCircuit.AllBusVmagPu; 
    VmaxMes(i) = max(AllBusVmagPu);

    totlss(i,:) = 1e-3*DSSObj.ActiveCircuit.Losses/sb; %pu. output seems to be in watts rather than kW
    totpwr(i,:) = DSSObj.ActiveCircuit.TotalPower/sb; %pu
end
toc

% withdraw results
DSSMon=DSSCircuit.Monitors; DSSMon.name='genvi';
DSSText.Command='show monitor genvi'; % this has to be in for some reason?

VmaxMat = reshape(VmaxMes,size(Pg));
VgenMat = reshape(VgenMes,size(Pg));

DSSMon=DSSCircuit.Monitors; DSSMon.name='genpq';
DSSText.Command='show monitor genpq';
Pgen_lin = -ExtractMonitorData(DSSMon,1,1)/sb; %pu, +ve is generation
Qgen_lin = -ExtractMonitorData(DSSMon,2,1)/sb; %pu, +ve is generation

%%
TotPwr = reshape(totpwr(:,1),size(Pg)) + 1i*reshape(totpwr(:,2),size(Pg)); %-ve implies load
TotLss = reshape(totlss(:,1),size(Pg)) + 1i*reshape(totlss(:,2),size(Pg)); %-ve implies load
Pgenmat=reshape(Pgen_lin,size(Pg));
Qgenmat=reshape(Qgen_lin,size(Pg));

dQ = 2e-2;
dV = 8e-3;

Qn = zeros(FF.n,numel(FF.Ps0_k));    e_L = zeros(FF.n,numel(FF.Ps0_k));
D_En = zeros(FF.n,numel(FF.Ps0_k));  D_Et = zeros(FF.n,numel(FF.Ps0_k));

[ Pb_fdr,Pp_fdr ] = find_pb_pp( Pgenmat,TotLss,VgenMat,Vp,1 );


tic
for i = 1:numel(FF.Ps0_k)
    i
    Ps = FF.Ps0*( Pb_fdr + (Pp_fdr - Pb_fdr)*FF.Ps0_k(i));
    [ Qn(:,i),D_En(:,i),D_Et(:,i),e_L(:,i) ] = calc_DE_EL( Ps,Pgenmat,Qgenmat,VmaxMat,VgenMat,Vp,FF.n,TotLss );
%     [ Qn(:,i),D_En(:,i),D_Et(:,i),e_L(:,i) ] = calc_DE_EL( FF.Ps0_k(i)*FF.Ps0,Pgenmat,Qgenmat,VmaxMat,VgenMat,Vp,FF.n,TotLss );
end
toc

tic
[ Qn_est,D_En_est,D_Et_est,e_L_est ] = estm_DE_EL( FF.Ps,FF.Ps0_k,Sload,Z,Vp,V0,20*FF.n );
toc

% subplot(211)
% plot(Qn/min(Qn),D_En); hold on;
% plot(Qn/min(Qn),D_Et,'--'); hold on;
% subplot(212)
% plot(Qn/min(Qn),D_Et./D_En); hold on;


% subplot(211)
% plot(Qn/min(Qn),D_En); hold on;
% plot(Qn/min(Qn),D_Et,'--'); hold on;
% subplot(212)
% plot(Qn/min(Qn),D_Et./D_En); hold on;


% %% calculate maximum power export s.t. voltage + current constraints:
% Vg_inV = VmaxMat + VNaN_outs; % remove numbers where the voltage is too high
% Vmn = (Vg_inV==max(Vg_inV));
% Vmax_pwr = max(real(TotPwr(Vmn)),[],2); % for each row, find the max power
% PgenV = Pgenmat(Vmn);
% 
% % Find 2 bus voltages, currents, reactive powers:
% [S0,Sg,Sl,~,Iest_pu,P0,Q0] = pred_S0_pscc( PgenV, Sload, Z, V0, Vp  );
% Iest = Iest_pu*ibN;
% 
% %% Return results RR
% [Pghat_est,Psub_hat_est] = lemma_1( Vp, Ip/ibN, V0, Z );
% [Pgprm_est,Psub_prm_est] = theorem_1( Vp, V0, Z );
% Pgen_hat_est = Pghat_est + real(Sload);
% Pgen_prm_est = Pgprm_est + real(Sload);
% 
% Psub_hat_meas = Imax_pwr(end);
% Pgen_hat_meas = Pgenmat(real(TotPwr)==Psub_hat_meas);
% Psub_prm_meas = max(Vmax_pwr);
% Pgen_prm_meas = Pgenmat(real(TotPwr)==Psub_prm_meas);
% 
% RR.Psub_hat = [Psub_hat_meas;Psub_hat_est];
% RR.Pgen_hat = [Pgen_hat_meas;Pgen_hat_est];
% RR.Psub_prm = [Psub_prm_meas;Psub_prm_est];
% RR.Pgen_prm = [Pgen_prm_meas;Pgen_prm_est];

RR.Qn = Qn;
RR.D_En = D_En;
RR.D_Et = D_Et;
RR.e_L = e_L;

RR.D_En_est = D_En_est;
RR.D_Et_est = D_Et_est;
RR.e_L_est = e_L_est;
RR.Qn_est = Qn_est;

RR.sbase = sb; RR.Z = Z; RR.Zabs = abs(Z); RR.lz = lambdas(Z); RR.V0 = V0;
RR.Ssc = Ssc; %kW 
RR.Sload = Sload; %in pu

%% Plotting options (pl_options)
figs = [];
% DAFS = 14; %default axis font size
% 
% if sum(ismember(pl_options,'pgp0'))  || sum(ismember(pl_options,'0'))
%     figs = [figs;figure('Color','White','Position',fig_nompos,'defaultaxesfontsize',DAFS,'defaulttextinterpreter','latex')];
%     plot(Pgenmat(Vmn),Vmax_pwr,'x-'); grid on; hold on;
%     plot(PgenV,P0);
%     axis equal;
%     axis([0 4 -1 1.5]);
%     xs = axis;
%     plot([1 1]*Pgen_prm_meas,xs(3:4),'k:');
% 
%     xlabel('$P_{gen}$ (pu)');
%     ylabel('$P_{0}^{Sub}$ (pu)');
%     lgnd=legend('Msrd.','Estd., (5)','$P_{g}''$, msrd.','Location','NorthWest');
%     set(lgnd,'interpreter','latex');
% 
% end
% if sum(ismember(pl_options,'pgqgq0')) || sum(ismember(pl_options,'0'))
%     figs = [figs;figure('Color','White','Position',fig_nompos,'defaultaxesfontsize',DAFS,'defaulttextinterpreter','latex')];
%     plot(Pgenmat(Vmn),Qgenmat(Vmn),'x-','Color',[0.3 0.3 0.3]); hold on; grid on;
%     plot(PgenV,imag(Sg + Sload),'k');
%     plot(Pgenmat(Vmn),imag(TotPwr(Vmn)),'x-','Color',[1.0 0.5 0.5]); hold on;
%     plot(PgenV,Q0,'r');
% 
%     lgnd = legend('$Q_{gen}$ msrd.','$Q_{gen}$ estd.','$Q_{comp}$ msrd.','$Q_{comp}$ estd.');
%     set(lgnd,'Interpreter','Latex');
%     xlabel('$P_{gen}$ (pu)'); ylabel('$Q$ (pu)'); axis equal;
% end
% if sum(ismember(pl_options,'imax')) || sum(ismember(pl_options,'0'))
%     figs = [figs;figure('Color','White','Position',fig_nompos,'defaultaxesfontsize',DAFS,'defaulttextinterpreter','latex')];%
%     plot(Pgenmat(Vmn),ImaxMat(Vmn),'x-'); grid on; hold on;
%     plot(PgenV,Iest); grid on; hold on;
%     
%     xs = axis;
%     plot(xs(1:2),[1 1]*FF.Ip,'k--');
%     plot([1 1]*Pgen_prm_meas,xs(3:4),'k:');
%     plot([1 1]*Pgen_hat_meas,xs(3:4),'k-.');
%     lgnd = legend('Max $I$, msrd.','$I_{g}$, estd.','$I_{+}$','$P_{g}''$, msrd.','$\hat{P}_{g}''$, msrd.','Location','NorthWest');
%     set(lgnd,'Interpreter','Latex');
%     xlabel('$P_{gen}$ (pu)'); ylabel('$|I|$ (A)');
% end

end

