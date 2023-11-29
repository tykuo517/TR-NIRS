%{
Check whether the fitting is complete

Benjamin Kao
Last update: 2021/01/18
%}

clc;clear;close all;

%% param
subject_name_arr={'KB'};
num_anser_to_generate=6; % number of target spec (true answer)
num_error_to_generate=1; % number of adding noise to the same, the first one will have no error
times_to_fitting=20; % number of fitting using different init value
input_dir='test_fitting_2023-04-17-09-33-13_copy'; % the folder containing the generated target spectrum
fitting_dir='fitting_SDS12345678910';

%% main
num_not_complete=0;

% for sbj=1:length(subject_name_arr)
for sbj=1
    for target_i=1:num_anser_to_generate
        for error_i=1:num_error_to_generate
            for i=1:times_to_fitting
                if exist(fullfile(input_dir,subject_name_arr{sbj},['target_' num2str(target_i) '_' num2str(error_i)],fitting_dir,['fitting_' num2str(i)],'fitted_spec.txt'),'file')==0
                    fprintf('%s target %d error %d fitting %d not exist!\n',subject_name_arr{sbj},target_i,error_i,i);
                    num_not_complete=num_not_complete+1;
                end
            end
        end
    end
end

fprintf('There are %d fitting not complete!\n',num_not_complete);