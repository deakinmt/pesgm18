close all; clear all; clc;
addpath('pesgm_funcs');

%%
% INPUTS TO DETERMINE BEHAVIOUR -------------------
C = {1.0,1.0,1.0,[1.0 0.6 0.2]};
VT = {1.04,1.08,1.12,1.05};

% WARNING! takes ~15 minutes per run (4x).
F.pg_ssc = linspace(-1e-6,0.2,1000);
F.qg_ssc = linspace(-0.4,1e-6,1200);
F.n = 120;
% F.pg_ssc = linspace(-1e-6,0.2,100);
% F.qg_ssc = linspace(-0.4,1e-6,100);
% F.n = 80;

%-------------------------------------------------

% remain the same:
F.SRC = 'SOURCEBUS';
F.filename = '\opendss_models\34Bus\ieee34Mod1_fxd';
F.feeder = '34bus'; % for calc_ww and T_ratios
F.tr_buses = {{'888','832'};{'832','852'};{'850','814'};{'800','SOURCEBUS'}};
F.fig_nompos = [200 150 450 340];

SCD = load('solar_curves_dataset.mat','Ps');
F.Ps0 = SCD.Ps; %pu
F.NUT = '834';
F.Vp = 1.06; %pu

RR = cell(size(C));
for i = 1:numel(C)
    F.V0 = VT{i};
    F.Ps0_k = C{i};
    RR{i} = run_pesgm_feeder( F );
end

% save('run_pesgm_34bus_rslt.mat');