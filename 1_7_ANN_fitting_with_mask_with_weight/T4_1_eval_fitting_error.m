%{
Arrange the fitting result of the random generated target spectrum

Ting-Yi Kuo
Last update: 2023/08/21
%}


clc;clear;close all; clearvars -global;

global lambda cw_net cw_param_range tr_net tr_param_range 

%% param
input_dir='test_fitting_2023-11-02-11-10-06'; % please move the fitting folders into this folder first.
subject_name_arr={'KB'}; %,'WH','WW'
num_anser_to_generate=15; % number of target spec (true answer)
num_error_to_generate=1; % number of adding noise to the same, the first one will have no error
times_to_fitting=20; % number of fitting using different init value

fitting_dir='fitting_SDS123456_234_gate1-5'; % the fitting folder

mode=1; % Fitting mode 1:TR+CW, 2:TR, 3:CW

% CW setting
num_SDS_cw=6; % how many SDS are in the target spectrum
SDS_choosed_cw=[1 2 3 4 5 6]; % the SDS chosen to fitted in the target spectrum
SDS_sim_correspond_cw=[1 2 3 4 5 6]; % the SDS index in the simulated spectrum corresponding to each SDS in the target spectrum, SDS_sim_correspond(x)=y means 'measure x = sim y'

% TR setting
num_SDS_tr=5;
SDS_choosed_tr=[2 3 4]; % the SDS chosen to fitted in the target spectrum
SDS_sim_correspond_tr=[1 2 3 4 5]; % the SDS index in the simulated spectrum corresponding to each SDS in the target spectrum, SDS_sim_correspond(x)=y means 'measure x = sim y'

num_gate=10;

% other setting
num_fitted_param=13; % the number of the fitted parameters

SDS_length_arr_cw=[0.8 1.5 2.12 3 3.35 4.5]; % cm
SDS_length_arr_tr=[1.5 2.2 2.9 3.6 4.3];

fitting_wl_file='fitting_wl.txt'; % the file in the epsilon folder containing the fitting wavelength
fitting_wl_tr=810;

model_dir='model_arrange'; % the folder containing the arranged ANN file

do_plot_anyPlot=1; % =1 to plot any plot
do_plot_individual_fitting=1; % =1 to plot the result of each individual fitting

%% init
fitting_wl=load(fullfile('epsilon',fitting_wl_file));
fitting_wl=[fitting_wl;fitting_wl_tr];
lambda=unique(fitting_wl);

mkdir(input_dir,'arrangement');

% load the param and OP answer
param_answer_arr=load(fullfile(input_dir,'answers','param_answer_arr.txt'));
OP_answer_arr=[];
for i=1:num_anser_to_generate
    temp_OP_arr=load(fullfile(input_dir,'answers',['OP_ans_' num2str(i) '.txt']));
    OP_answer_arr(:,:,i)=interp1(temp_OP_arr(:,1),temp_OP_arr(:,2:end),lambda);
end
% OP_answer_arr=OP_answer_arr(:,[3 5 9],:);
% OP_answer_arr=OP_answer_arr(:,[2 3 4],:);

%% main

for sbj=1:length(subject_name_arr)
    % make output folder and find fitting folder
%     mkdir(fullfile(input_dir,'arrangement',subject_name_arr{sbj}));
    
    % load ANN model
    load(fullfile(model_dir,[subject_name_arr{sbj} '_cw_model.mat'])); % cw_net, cw_param_range
    load(fullfile(model_dir,[subject_name_arr{sbj} '_tr_model.mat'])); % tr_net, tr_param_range
    
    index=0;
    for target_i=1:num_anser_to_generate
        for error_i=1:num_error_to_generate
            index=index+1;
%             fprintf('Process %s target %d, error %d, fitting ',subject_name_arr{sbj},target_i,error_i);

            temp_target_folder=fullfile(input_dir,subject_name_arr{sbj},['target_' num2str(target_i) '_' num2str(error_i)],fitting_dir);

            target_spec=load(fullfile(input_dir,subject_name_arr{sbj},['target_' num2str(target_i) '_' num2str(error_i)],'target_spec.txt'));
            target_dtof=load(fullfile(input_dir,subject_name_arr{sbj},['target_' num2str(target_i) '_' num2str(error_i)],'target_dtof.txt'));

            fitting_result_sort=load(fullfile(temp_target_folder,'fitting_result_sort.txt'));
            best_result_index=fitting_result_sort(1,1);
            load(fullfile(temp_target_folder,'mask.mat'));
            weight=load(fullfile(temp_target_folder,'weight.txt'));
            
            fitting_op_arr(:,:,index)=load(fullfile(temp_target_folder,['fitting_' num2str(best_result_index)],'fitted_mu.txt'));
            fitted_ann_spec_arr(:,:,index)=load(fullfile(temp_target_folder,['fitting_' num2str(best_result_index)],'fitted_spec.txt'));
            fitted_ann_dtof_arr_2D(:,:,index)=load(fullfile(temp_target_folder,['fitting_' num2str(best_result_index)],'fitted_dtof.txt'));

            fun_init_param_to_mu_spec(); % load the epsilon for fitting wl

            interped_target_spec=interp1(target_spec(:,1),target_spec(:,2:end),lambda); % for calculate error
 
            % calculate the error
            % for CW
            calculated_fitted_spec_error(:,:,index)=sqrt(mean((fitted_ann_spec_arr(:,:,index)./interped_target_spec-1).^2,1));

            % for TR
            calculated_fitted_error_tr=fitted_ann_dtof_arr_2D(:,:,index)./target_dtof(:,:)-1;
            calculated_fitted_dtof_error(:,:,index)=calculated_fitted_error_tr; % save errors
        end
    end
    
    % calculate mean error of each SDS and gate
    mean_fitted_dtof_error=mean(calculated_fitted_dtof_error,3);
    
    select_gate=sum(mask,1);
end
