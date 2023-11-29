%{
Arrange the training error of each subject and each SDS

Benjamin Kao
Last update: 2020/01/04
%}

clc;clear;close all;

%% param
subject_name_arr={'ZJ','YH','YF','WW','WH','SJ','SC','KB','BT'}; % the name of the subject
input_dir_arr={'ZJ_2020-12-22-12-36-51','YH_2020-12-22-22-44-50','YF_2020-12-22-15-14-33','WW_2020-12-23-18-41-36','WH_2020-12-22-20-37-29','SJ_2020-12-22-18-48-01','SC_2020-12-23-02-22-28','KB_2020-12-23-16-36-17','BT_2020-12-23-00-35-30'};
num_SDS=7;

%% main

total_rmspe_arr=[];
for i=1:length(input_dir_arr)
    temp_error=load(fullfile(input_dir_arr{i},'error_record.txt'));
    total_rmspe_arr(end+1,:)=temp_error(:,3)';
end

for i=1:size(total_rmspe_arr,1)
    for j=1:size(total_rmspe_arr,2)
        fprintf('\t%.2f%%',total_rmspe_arr(i,j)*100);
    end
    fprintf('\n');
end