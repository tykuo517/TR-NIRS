%{
Generate the database for the fitting program to find the initial value

Benjamin Kao
Last update: 2021/01/19
%}

clc;clear;close all;clearvars -global;

global lambda fitting_wl_tr Lbound Ubound cw_net cw_param_range tr_net tr_param_range;

%% param
model_dir='model_arrange'; % the folder of the models
subject_name_arr={'KB'}; % the name of the subjects
num_SDS_cw=7;
num_SDS_tr=5;
num_gate=10; % the number of SDS, in the ANN

num_init_spec=200000; % number of random generated spectrum
fitting_wl=(700:4:900)'; % the wavelenth to generate the spectrum
fitting_wl_tr=810;

fitting_wl_file='fitting_wl.txt'; % the file in the epsilon folder containing the fitting wavelength

output_dir='initValue_database_2'; % the folder to save the output database

% the range for the parameters
%                      hc_1    sto2_1  hc_2    sto2_2 hc_4     sto2_4  mel
mua_param_Lbound=     [1       0.3     1       0.3    50       0.3     0];
mua_param_Ubound=     [150     1       150     1      200      1       0.003 ];

% %                      hc_2    sto2_2  hc_4  sto2_4 
% mua_param_Lbound=     [1       0.3     50      0.3];
% mua_param_Ubound=     [150     1       200     1 ];


%% init
mkdir(output_dir);

rng('shuffle'); % random seed

wl_used_to_fitting=load(fullfile('epsilon',fitting_wl_file));

% find all the wl to check
wl_to_check=[wl_used_to_fitting; fitting_wl];
wl_to_check=unique(wl_to_check);

%% main

for sbj_i=1:length(subject_name_arr)
    lambda=wl_to_check;
    fun_init_param_to_mu_spec(); % load the epsilon for target wl
    
    % load ANN model
    load(fullfile(model_dir,[subject_name_arr{sbj_i} '_cw_model.mat'])); % cw_net, cw_param_range
    load(fullfile(model_dir,[subject_name_arr{sbj_i} '_tr_model.mat'])); % tr_net, tr_param_range
    
    % load AK range
    A_Krange=load(fullfile(model_dir,['A_Krange_arr_' subject_name_arr{sbj_i} '.mat']));
    Lbound=[A_Krange.Lbound mua_param_Lbound]; % add the bound for A, K for L1, 2, 4
    Ubound=[A_Krange.Ubound mua_param_Ubound];


    % generate the init value array
    init_value_arr=zeros(num_init_spec,size(Lbound,2));
    tic;
    for init_i=1:num_init_spec

        fprintf('Finding init value %d, ',init_i);

        %% generate init value and spec 先根據Upper,Lower bound 隨機取值,再將散射係數的值找到K值
        temp_init_value=rand(size(Lbound)).*(Ubound-Lbound)+Lbound; % generate the init param
        % adjust the K range by the A value
        for L=1:3
            temp_Krange=interp1(A_Krange.A_Krange_arr{L}(:,1),A_Krange.A_Krange_arr{L}(:,2:3),temp_init_value(2*L-1),'pchip');
            temp_init_value(2*L)=rand(1,1).*(temp_Krange(1)-temp_Krange(2))+temp_Krange(2);
        end
        [OP_arr,~]=fun_param_to_mu(temp_init_value,0);

        retry_count=1;
        while ~fun_in_OP_range(OP_arr)
            fprintf('%d ',retry_count); retry_count=retry_count+1;
            temp_init_value=rand(size(Lbound)).*(Ubound-Lbound)+Lbound;
            for L=1:3
                temp_Krange=interp1(A_Krange.A_Krange_arr{L}(:,1),A_Krange.A_Krange_arr{L}(:,2:3),temp_init_value(2*L-1),'pchip');
                temp_init_value(2*L)=rand(1,1).*(temp_Krange(1)-temp_Krange(2))+temp_Krange(2);
            end
            [OP_arr,~]=fun_param_to_mu(temp_init_value,0);
        end

        init_value_arr(init_i,:)=temp_init_value;
        fprintf('\n');
    end
    procsee_time=toc;

    % forward the spectrum
    lambda=fitting_wl;
    fun_init_param_to_mu_spec(); % load the epsilon for target wl
    
    init_spec_arr=zeros(length(fitting_wl),num_SDS_cw,num_init_spec);
    init_dtof_arr=zeros(num_gate,num_SDS_tr,num_init_spec);

    % for each target spec
    for init_i=1:num_init_spec
        fprintf('Subject %d spec %d\n',sbj_i,init_i);

        %% generate init value and spec
        [OP_arr,~]=fun_param_to_mu(init_value_arr(init_i,:),0);

        init_spec_arr(:,:,init_i)=fun_ANN_forward(OP_arr,0);
        
        temp_OP_arr=[fitting_wl OP_arr];
        OP_arr=interp1(temp_OP_arr(:,1),temp_OP_arr(:,2:end),fitting_wl_tr,'pchip');
        
        target_dtof_=fun_ANN_forward(OP_arr,1);
        
        for i=1:num_SDS_tr
            target_dtof(:,i)=target_dtof_((i-1)*num_gate+1:i*num_gate)';
        end
        
        init_dtof_arr(:,:,init_i)=target_dtof;
        
    end
    
    save(fullfile(output_dir,[subject_name_arr{sbj_i} '_DB.mat']),'init_spec_arr','init_dtof_arr','init_value_arr','fitting_wl','fitting_wl_tr');
end

disp('Done!');