%{
Find the noise add to the target spctrum

Benjamin Kao
Last update: 2021/01/19
%}

clc;clear;close all; clearvars -global;

%% param
input_dir='test_fitting_2023-10-12-00-56-05'; % please move the fitting folders into this folder first.
subject_name_arr={'KB'};
num_anser_to_generate=15; % number of target spec (true answer)
num_error_to_generate=15; % number of adding noise to the same, the first one will have no error

%% main
for sbj=1:length(subject_name_arr)
    SDS_added_noise_arr=[];
    for target_i=1:num_anser_to_generate
        fprintf('Processing %s target %d\n',subject_name_arr{sbj},target_i);
        for error_i=1:num_error_to_generate
            temp_spec=load(fullfile(input_dir,subject_name_arr{sbj},['target_' num2str(target_i) '_' num2str(error_i)],'target_spec.txt'));
            temp_spec=temp_spec(:,2:end); % delete the wl column
%             temp_spec=temp_spec(101,2:end); % use the 800nm reflectance
            if error_i==1
                baseline_value=temp_spec;
            end
            SDS_added_noise_arr(error_i,:,target_i)=sqrt(mean((temp_spec./baseline_value-1).^2,1));
        end
    end
    save(fullfile(input_dir,subject_name_arr{sbj},'SDS_added_noise_arr.mat'),'SDS_added_noise_arr');
end