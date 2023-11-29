%{
Generate the target param answer, to let the all subject have the same answer to fitting
Also add the noise to the simultion spectrum and DTOF

Ting-Yi Kuo
Last update: 2023/08/21
%}

clc;clear;close all;

global lambda fitting_wl_tr Lbound Ubound cw_net tr_net cw_param_range tr_param_range;

%% param
subject_name_arr={'KB'};%,'WH','WW'
num_anser_to_generate=15; % number of target spec (true answer)
num_error_to_generate=15; % number of adding noise to the same, the first one will have no error

num_SDS_cw=6;
num_SDS_tr=5;
num_gate=10;

% SDS_error_arr=[2 2.3 0.75 4 1 2]; % in %, the error (CV) of each SDS
SDS_error_arr=[3 4.2 5.1 5.2 5.4 12.1]; % in %, the error (CV) of each SDS

% load('CV_10_11.mat');
load('../1_3_MCX_lookup_table/results/cv_1E11'); % cv_value
gate_error_arr=100*cv_value(1:num_gate,1:num_SDS_tr);

use_slope_mode=1; % if =1, use slope noise instead of constant noise
small_noise_ratio=0.05;


target_spec_wl=(700:900)';
fitting_wl_tr=810;

%% init
%            hc_1    sto2_1  hc_2    sto2_2 hc_4     sto2_4  mel
Lbound=     [1       0.3     1       0.3    50       0.3     0];
Ubound=     [150     1       150     1      200      1       0.003 ];
% %            hc_2    sto2_2 hc_4     sto2_4 
% Lbound=     [1       0.3    50       0.3 ];
% Ubound=     [150     1      200      1    ];

model_dir='model_arrange';
A_Krange=load(fullfile(model_dir,'A_Krange_arr_KB.mat'));
Lbound=[A_Krange.Lbound Lbound]; % add the bound for A, K for L1, 2, 4
Ubound=[A_Krange.Ubound Ubound];

rng('shuffle'); % random seed
fitting_wl=[target_spec_wl;fitting_wl_tr];
lambda=unique(fitting_wl);

fun_init_param_to_mu_spec(); % load the epsilon for target wl

output_dir=['test_fitting_' datestr(datetime('now'),'yyyy-mm-dd-HH-MM-ss')];
mkdir(output_dir);
mkdir(fullfile(output_dir,'answers'));

%% main

%% generate the target answer using sbj 1 (small mus range)

sbj=1;
% load ANN model, for the param range
load(fullfile(model_dir,[subject_name_arr{sbj} '_cw_model.mat'])); % cw_net, cw_param_range
load(fullfile(model_dir,[subject_name_arr{sbj} '_tr_model.mat'])); % tr_net, tr_param_range

% generate answers
param_answer_arr=[];
for target_i=1:num_anser_to_generate
    fprintf('Generating answer %d\n',target_i);
    temp_init_value=rand(size(Lbound)).*(Ubound-Lbound)+Lbound;
    for L=1:3
        temp_Krange=interp1(A_Krange.A_Krange_arr{L}(:,1),A_Krange.A_Krange_arr{L}(:,2:3),temp_init_value(2*L-1),'pchip');
        temp_init_value(2*L)=rand(1,1).*(temp_Krange(1)-temp_Krange(2))+temp_Krange(2);
    end
    [OP_arr,~]=fun_param_to_mu(temp_init_value,0);

    while ~fun_in_OP_range(OP_arr)
        temp_init_value=rand(size(Lbound)).*(Ubound-Lbound)+Lbound;
        for L=1:3
            temp_Krange=interp1(A_Krange.A_Krange_arr{L}(:,1),A_Krange.A_Krange_arr{L}(:,2:3),temp_init_value(2*L-1),'pchip');
            temp_init_value(2*L)=rand(1,1).*(temp_Krange(1)-temp_Krange(2))+temp_Krange(2);
        end
        [OP_arr,~]=fun_param_to_mu(temp_init_value,0);
    end

    param_answer_arr(target_i,:)=temp_init_value;
end
save(fullfile(output_dir,'answers','param_answer_arr.txt'),'param_answer_arr','-ascii','-tabs');

%% generate target spec
fprintf('Generating target spec.\n');
for sbj=1:length(subject_name_arr)
    fprintf('\tSubject %s\n',subject_name_arr{sbj});
    
    for target_i=1:num_anser_to_generate
        % generate OP
        [OP_arr,~]=fun_param_to_mu(param_answer_arr(target_i,:),0);
        
        if sbj==1
            to_save=[target_spec_wl OP_arr];
            save(fullfile(output_dir,'answers',['OP_ans_' num2str(target_i) '.txt']),'to_save','-ascii','-tabs');
        end
        
        % generate original spectral
        orig_target_spec=fun_ANN_forward(OP_arr,0);
        
        OP_arr=interp1(target_spec_wl,OP_arr,fitting_wl_tr,'pchip');
        orig_target_dtof_=fun_ANN_forward(OP_arr,1);
        
        for i=1:num_SDS_tr
            orig_target_dtof(:,i)=orig_target_dtof_((i-1)*num_gate+1:i*num_gate)';
        end
        
        % add noise
        for error_i=1:num_error_to_generate
            mkdir(fullfile(output_dir,subject_name_arr{sbj},['target_' num2str(target_i) '_' num2str(error_i)]));
            % for CW
            if error_i==1
                temp_SDS_error=zeros(1,num_SDS_cw);
            else
                if use_slope_mode==1
                    temp_SDS_error=[];
                    for s=1:num_SDS_cw
                        % decide the 2 end of the slope line
                        error_end=normrnd(0,SDS_error_arr(s),1,2);
                        error_end(2)=error_end(1)+(error_end(2)-error_end(1))*rand(1);
                        SDS_error_spec=interp1([min(target_spec_wl) max(target_spec_wl)],error_end,target_spec_wl);

                        % add small noise
                        SDS_error_spec=SDS_error_spec+normrnd(0,SDS_error_arr(s)*small_noise_ratio,size(target_spec_wl));
                        temp_SDS_error(:,s)=SDS_error_spec;
                    end
                else
                    temp_SDS_error=normrnd(0,SDS_error_arr);
                end
            end
            temp_target_spec=orig_target_spec(:,1:num_SDS_cw).*(1+temp_SDS_error/100);
            to_save=[target_spec_wl temp_target_spec];
            save(fullfile(output_dir,subject_name_arr{sbj},['target_' num2str(target_i) '_' num2str(error_i)],'target_spec.txt'),'to_save','-ascii','-tabs');
            
            % for TR
            if error_i==1
                temp_SDS_error=zeros(1,num_SDS_tr);
            else
                temp_SDS_error=normrnd(0,gate_error_arr);
            end
            temp_target_dtof=orig_target_dtof(:,1:num_SDS_tr).*(1+temp_SDS_error/100);
            to_save=[temp_target_dtof];
            save(fullfile(output_dir,subject_name_arr{sbj},['target_' num2str(target_i) '_' num2str(error_i)],'target_dtof.txt'),'to_save','-ascii','-tabs');
        end
    end
end

disp('Done!');