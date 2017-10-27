close all; clear all; clc;
cd('C:\Users\chri3793\Documents\MATLAB\DPhil\pesgm18_mtlb')
% cd('C:\Users\Matt\Documents\MATLAB\DPhil\pesgm18_mtlb');
addpath('pesgm_funcs');
% fig_loc = [pwd,'\figures\'];
% fig_loc = 'C:\Users\chri3793\Documents\MATLAB\DPhil\pesgm18_mtlb\figures\';
fig_loc = 'C:\Users\Matt\Documents\MATLAB\DPhil\pesgm18_mtlb\figures\';
%%
cd('C:\Users\chri3793\Documents\MATLAB\DPhil\pesgm18_mtlb')
% cd('C:\Users\Matt\Documents\MATLAB\DPhil\pesgm18_mtlb')
% INPUTS TO DETERMINE BEHAVIOUR -------------------
F.Vp = 1.06; %pu
% F.pg_ssc = linspace(-1e-3,0.2,50);
% F.qg_ssc = linspace(-0.4,1e-3,50);

% F.pg_ssc = linspace(-1e-3,0.15,100);
% F.qg_ssc = linspace(-0.1,1e-3,100);

F.pg_ssc = linspace(-1e-3,0.2,100);
F.qg_ssc = linspace(-0.4,1e-3,300);

% remain the same:
F.SRC = 'SOURCEBUS';

F.filename = '\opendss_models\34Bus\ieee34Mod1_fxd';
F.feeder = '34bus'; % for calc_ww and T_ratios
F.tr_buses = {{'888','832'};{'832','852'};{'850','814'};{'800','SOURCEBUS'}};
F.fig_nompos = [200 150 450 340];

SCD = load('solar_curves_dataset.mat','Ps');
F.Ps = 2.7*SCD.Ps; %pu
F.Qmaxpu = 4;

% First results figure: 834 bus feeder voltage
figname1 = [fig_loc,'34bus_pgp0'];
F.NUT = '834';

pl_options={'pgp0'};
[ ~,figs ] = run_pesgm_feeder( F,pl_options );