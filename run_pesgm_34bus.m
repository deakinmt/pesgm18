close all; clear all; clc;
cd('C:\Users\chri3793\Documents\MATLAB\DPhil\pesgm18_mtlb')
% cd('C:\Users\Matt\Documents\MATLAB\DPhil\pscc18_mtlb');
addpath('pesgm_funcs');
% fig_loc = [pwd,'\figures\'];
fig_loc = 'C:\Users\chri3793\Documents\MATLAB\DPhil\pesgm18_mtlb\figures\';
% fig_loc = 'C:\Users\Matt\Documents\DPhil\pscc18\pscc18_paper\figures\';
%%
cd('C:\Users\chri3793\Documents\MATLAB\DPhil\pesgm18_mtlb')
% INPUTS TO DETERMINE BEHAVIOUR -------------------
F.Vp = 1.06; %pu
% F.pg_ssc = linspace(-0.001,0.2,20);
% F.qg_ssc = linspace(-0.4,0.03,20);
F.pg_ssc = linspace(-0.001,0.2,100);
F.qg_ssc = linspace(-0.4,0.03,300);

% remain the same:
F.SRC = 'SOURCEBUS';

F.filename = '\opendss_models\34Bus\ieee34Mod1_fxd';
F.feeder = '34bus'; % for calc_ww and T_ratios
F.tr_buses = {{'888','832'};{'832','852'};{'850','814'};{'800','SOURCEBUS'}};
F.fig_nompos = [200 150 450 340];

% First results figure: 834 bus feeder voltage
figname1 = [fig_loc,'34bus_pgp0'];
F.NUT = '834';
F.pg_ssc = linspace(0.001,0.2,20);
F.qg_ssc = linspace(-0.4,0.03,20);
pl_options={'pgp0'};
[ ~,figs ] = run_pesgm_feeder( F,pl_options );