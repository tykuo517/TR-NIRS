%{
copy the target spec (without fitting dirs) to make new fitting

Benjamin Kao
Last update: 2021/03/17
%}

clc;clear;close all;

global lambda Lbound Ubound net param_range;

%% param
subject_name_arr={'KB'};
num_anser_to_generate=6; % number of target spec (true answer)
num_error_to_generate=1; % number of adding noise to the same, the first one will have no error
input_dir='test_fitting_2023-04-10-21-23-21';
output_dir='test_fitting_2023-04-10-21-23-21_copy';

%% main
mkdir(output_dir);

% for sbj=1:length(subject_name_arr)
for sbj=1
    fprintf('sbj %d\n',sbj);
    for i=1:num_anser_to_generate
        for j=1:num_error_to_generate
            mkdir(fullfile(output_dir,subject_name_arr{sbj},['target_' num2str(i) '_' num2str(j)]));
            copyfile(fullfile(input_dir,subject_name_arr{sbj},['target_' num2str(i) '_' num2str(j)],'target_spec.txt'),fullfile(output_dir,subject_name_arr{sbj},['target_' num2str(i) '_' num2str(j)],'target_spec.txt'));
        end
    end
end

disp('Done!');