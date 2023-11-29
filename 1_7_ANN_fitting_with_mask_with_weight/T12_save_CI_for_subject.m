%{
Save the confidence interval for each subject

Benjamin Kao
Last update: 2021/03/17
%}

clc;clear;close all; clearvars -global;

%% param
output_dir='fitted_result';
fitting_dir_arr={'test_fitting_2021-01-17-16-47-46','test_fitting_2021-01-22-03-42-40'}; % the mother dir of the testing target fitting
fitting_SDS_dir_arr={'fitting_SDS1234','fitting_SDS1234','fitting_SDS346','fitting_SDS2345','fitting_SDS1234','fitting_SDS12345','fitting_SDS123456'}; % the sub fitting dir in the testing target fitting dir for each SDS combination
fitting_dir_index=[1 1 2 2 2 2 2];
add_noise_arr=[0 1 1 1 1 1 1];
fitting_name_arr={'SDS 1234 merged no noise','SDS 1234 add system noise','SDS 346 error','SDS 2345 error','SDS 1234 error','SDS 12345 error','SDS 123456 error'}; % the name for each SDS or error condition, just for understandin
target_name_arr={'tc_mean','ww_mean','wh_mean','yf_mean','kb_mean'}; % the name of the target spectrum
subject_fitting_SDS_arr=[7 3 5 6 5];

%% main
for sbj=1:length(target_name_arr)
    i=subject_fitting_SDS_arr(sbj);
    if add_noise_arr(i)==0
        to_load_name='OP_error_arr_noError.mat';
    else
        to_load_name='OP_error_arr_Error.mat';
    end
    error_info=load(fullfile(fitting_dir_arr{fitting_dir_index(i)},'arrangement',fitting_SDS_dir_arr{i},to_load_name));
    error_CI=error_info.OP_error_CI;
    save(fullfile(output_dir,[target_name_arr{sbj} '_error_CI.mat']),'error_CI');
end