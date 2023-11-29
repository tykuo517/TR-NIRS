clc;clear;close all;
%% init
CS_param=load('CS_baseline_param.txt');

global lambda;
lambda=CS_param(:,1);

param = containers.Map;

% epidermis with melanin, from skin P12 result
param('mel_1') = 0.016502;
param('A1_1') = 4584.67793;
param('K1_1') = 2.186113;
param('A2_1') = 499.912945;
param('K2_1') = 0.370345;
param('th_1') = 0.010942;
param('g_1') = 0.75;
% dermis, from skin P12 result
param('A1_2') = 19510.9413;
param('K1_2') = 2.998319;
param('A2_2') = 288.304257;
param('K2_2') = 0.364456;
param('hc_2') = 0.001502;
param('sto2_2') = 0.690279;
param('th_2') = 0.6 - param('th_1');
param('g_2') = 0.715;
% skull, from tin CS result
param('sto2_3') = 0.58964876;
param('hc_3') = 54.797124;
param('A_3') = 1.3901313;
param('K_3') = 0.98902308;
param('th_3') = 0.75;
param('g_3') = 0.9;
% CSF
param('th_4') = 0.6;
param('g_4') = 0.9;
% gray matter, from tin CS result
param('A_5') = 1507230.5 * 10; % turn 1/mm into 1/cm
param('K_5') = 3.0906513;
param('sto2_5') = 0.67456;
param('hc_5') = 92.7674;
param('th_5') = 0.8;
param('g_5') = 0.9;

%% fitting layer 2
x0=[param('hc_2') param('sto2_2')];
new_x=fminsearch();

%% function
function error=layer2mua(x0)
    mu=param_to_mu_spec(param);
end