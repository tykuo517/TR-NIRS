%{
Fitting the target spec using ANN

Benjamin Kao
Last update: 2020/01/25
%}

clc;clear;close all;

global lambda Lbound Ubound net param_range target_spec orig_target_spec num_SDS SDS_choosed A_Krange SDS_sim_correspond;

%% param
target_name_arr={'tc_mean','ww_mean','wh_mean','yf_mean','kb_mean'}; % the name of the target spectrum
model_name_arr={'ZJ','WW2','WH2','YF','KB'}; % the name of ANN model corresponding to each target spec

output_folder='fitting_SDS346_3'; % the name of output main folder

num_SDS=6; % how many SDS are in the target spectrum
SDS_choosed=[3 4 6]; % the SDS chosen to fitted in the target spectrum
SDS_sim_correspond=[1 2 3 4 5 6]; % the SDS index in the simulated spectrum corresponding to each SDS in the target spectrum, SDS_sim_correspond(x)=y means 'measure x = sim y'

times_to_fitting=20; % number of init value used to fitting

fitting_wl_file='fitting_wl.txt'; % the file in the epsilon folder containing the fitting wavelength

model_dir='model_arrange'; % the folder containing the arranged ANN file
target_dir='input_target_2'; % the folder containing the target spectrum

%% init
assert(num_SDS==length(SDS_sim_correspond),'Error: ''SDS_sim_correspond'' not match to ''num_SDS''!');

%                      hc_1    sto2_1 hc_2    sto2_2 hc_4     sto2_4   mel
mua_param_Lbound=     [1       0.3     1       0.3    50       0.3     0];
mua_param_Ubound=     [150     1       150     1      200      1       0.003 ];

mua_coef_bound=[mua_param_Lbound; mua_param_Ubound];

fitting_wl=load(fullfile('epsilon',fitting_wl_file));
lambda=fitting_wl;
fun_init_param_to_mu_spec(); % load the epsilon for fitting wl

mkdir(output_folder);
save(fullfile(output_folder,'mua_coef_bound.txt'),'mua_coef_bound','-ascii','-tabs');

%% main

for sbj=1:length(target_name_arr)
    A_Krange=load(fullfile(model_dir,['A_Krange_arr_' model_name_arr{sbj} '.mat']));
    Lbound=[A_Krange.Lbound mua_param_Lbound]; % add the bound for A, K for L1, 2, 4
    Ubound=[A_Krange.Ubound mua_param_Ubound];
    
    output_dir=fullfile(output_folder,[target_name_arr{sbj} '_' datestr(datetime('now'),'yyyy-mm-dd-HH-MM-ss')]);
    mkdir(output_dir);
    
    % backup this script
    copyfile([mfilename '.m'],fullfile(output_dir,[mfilename '.m']));
    
    % load target spec, and interpolation
    target_spec=load(fullfile(target_dir,[target_name_arr{sbj} '.txt']));
    orig_target_spec=target_spec;
    target_spec=interp1(target_spec(:,1),target_spec(:,2:end),fitting_wl);
    
    % load ANN model
    ANN_model=load(fullfile(model_dir,[model_name_arr{sbj} '_model.mat'])); % net, param_range
    net=ANN_model.net; param_range=ANN_model.param_range;
        
    %% fitting
    % find the init value
    [init_arr,init_error]=fun_choose_initValue(orig_target_spec(:,2:end),orig_target_spec(:,1),model_name_arr{sbj});
    init_arr=init_arr(1:times_to_fitting,:);

    init_error=init_error(1:times_to_fitting);
    to_save=[init_arr init_error];
    
    fid=fopen(fullfile(output_dir,'init_arr.csv'),'a');
    fprintf(fid,'A1,K1,A2,K2,A4,K4,hc_1,StO2_1,hc_2,StO2_2,hc_4,StO2_4,mel,error,\n');
    for fitting_i=1:times_to_fitting
        fprintf(fid,'%f,',to_save(fitting_i,:));
        fprintf(fid,'\n');
    end
    fclose(fid);
    save(fullfile(output_dir,'init_arr.txt'),'to_save','-ascii','-tabs');
    
    % start fitting using different init value
    for fitting_i=1:times_to_fitting
        mkdir(fullfile(output_dir,['fitting_' num2str(fitting_i)]));

        fun_auto_fittingMuaMus(init_arr(fitting_i,:),fullfile(output_dir,['fitting_' num2str(fitting_i)]));
    end
        
end

%% function