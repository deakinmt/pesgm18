function RR = run_pesgm_feeder( FF )

Vp = FF.Vp;

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

Qgd = zeros(FF.n,numel(FF.Ps0_k));   e_L = zeros(FF.n,numel(FF.Ps0_k));
D_Eg = zeros(FF.n,numel(FF.Ps0_k));  D_Et = zeros(FF.n,numel(FF.Ps0_k));

[ Pb_fdr,Pp_fdr ] = find_pb_pp( Pgenmat,TotLss,VgenMat,Vp );

tic
for i = 1:numel(FF.Ps0_k)
    i
    Ps = FF.Ps0*( Pb_fdr + (Pp_fdr - Pb_fdr)*FF.Ps0_k(i));
    [ Qgd(:,i),D_Eg(:,i),D_Et(:,i),e_L(:,i) ] = calc_DE_EL( Ps,Pgenmat,Qgenmat,VmaxMat,VgenMat,Vp,FF.n,TotLss );
end
toc

tic
[ Qgd_est,D_Eg_est,D_Et_est,e_L_est ] = estm_DE_EL( FF.Ps0,FF.Ps0_k,Sload,Z,Vp,V0,FF.n );
toc

RR.Qgd = Qgd;
RR.D_Eg = D_Eg;
RR.D_Et = D_Et;
RR.e_L = e_L;

RR.Qgd_est = Qgd_est;
RR.D_Eg_est = D_Eg_est;
RR.D_Et_est = D_Et_est;
RR.e_L_est = e_L_est;


RR.sbase = sb; RR.Z = Z; RR.Zabs = abs(Z); RR.lz = lambdas(Z); RR.V0 = V0;
RR.Ssc = Ssc; %kW 
RR.Sload = Sload; %in pu

end

