%{
Process the multiple solution subjects, choose the solution(s) to save

Benjamin Kao
Last update: 2021/03/10
%}

clc;clear;close all; clearvars -global;

global lambda net param_range 

%% param
input_dir='fitted_result';
target_name_arr={'tc_mean','ww_mean','wh_mean','yf_mean','kb_mean'}; % the name of the target spectrum
to_process_sbj_index=3; % the index of the subject to process
to_save_multiSol_index=[1]; % the array containing the index of the solution(s) to save

multiSol_dir='multiSol';

%% init 
if exist(fullfile(input_dir,multiSol_dir),'dir')==0
    mkdir(fullfile(input_dir,multiSol_dir));
end

%% main

% copy the old files
sbj_name=target_name_arr{to_process_sbj_index};
choosed_rank=load(fullfile(input_dir,[sbj_name '_choosed_rank.txt']));

assert(length(choosed_rank)>=max(to_save_multiSol_index));

for i=1:length(choosed_rank)
    try
        movefile(fullfile(input_dir,[sbj_name '_fitted_OP_' num2str(i) '.txt']),fullfile(input_dir,multiSol_dir,[sbj_name '_fitted_OP_' num2str(i) '.txt']));
    catch
        fprintf('%s don''t exist!\n',[sbj_name '_fitted_OP_' num2str(i) '.txt']);
    end
    try
        movefile(fullfile(input_dir,[sbj_name '_fitted_spec_' num2str(i) '.txt']),fullfile(input_dir,multiSol_dir,[sbj_name '_fitted_spec_' num2str(i) '.txt']));
    catch
        fprintf('%s don''t exist!\n',[sbj_name '_fitted_spec_' num2str(i) '.txt']);
    end
end

try
    movefile(fullfile(input_dir,[sbj_name '_fitted_OP_info.mat']),fullfile(input_dir,multiSol_dir,[sbj_name '_fitted_OP_info.mat']));
catch
    fprintf('%s don''t exist!\n',[sbj_name '_fitted_OP_info.mat']);
end

try
    movefile(fullfile(input_dir,[sbj_name '_choosed_rank.txt']),fullfile(input_dir,multiSol_dir,[sbj_name '_choosed_rank.txt']));
catch
    fprintf('%s don''t exist!\n',[sbj_name '_choosed_rank.txt']);
end

try
    movefile(fullfile(input_dir,[sbj_name '_choosed_spec.png']),fullfile(input_dir,multiSol_dir,[sbj_name '_choosed_spec.png']));
catch
    fprintf('%s don''t exist!\n',[sbj_name '_choosed_spec.png']);
end

try
    movefile(fullfile(input_dir,[sbj_name '_choose_process.png']),fullfile(input_dir,multiSol_dir,[sbj_name '_choose_process.png']));
catch
    fprintf('%s don''t exist!\n',[sbj_name '_choose_process.png']);
end

try
    movefile(fullfile(input_dir,[sbj_name '_choosed_OP.png']),fullfile(input_dir,multiSol_dir,[sbj_name '_choosed_OP.png']));
catch
    fprintf('%s don''t exist!\n',[sbj_name '_choosed_OP.png']);
end

% save the new files
new_choosed_rank=choosed_rank(to_save_multiSol_index);
save(fullfile(input_dir,[sbj_name '_choosed_rank.txt']),'new_choosed_rank','-ascii','-tabs');

load(fullfile(input_dir,multiSol_dir,[sbj_name '_fitted_OP_info.mat']));
fitted_OP_arr=fitted_OP_arr(:,:,to_save_multiSol_index);
toOutput_rank=new_choosed_rank;
save(fullfile(input_dir,[sbj_name '_fitted_OP_info.mat']),'OP_CV','fitted_OP_arr','fitting_index','lambda','toOutput_rank');

for i=1:length(to_save_multiSol_index)
    copyfile(fullfile(input_dir,multiSol_dir,[sbj_name '_fitted_OP_' num2str(to_save_multiSol_index(i)) '.txt']),fullfile(input_dir,[sbj_name '_fitted_OP_' num2str(i) '.txt']));
    copyfile(fullfile(input_dir,multiSol_dir,[sbj_name '_fitted_spec_' num2str(to_save_multiSol_index(i)) '.txt']),fullfile(input_dir,[sbj_name '_fitted_spec_' num2str(i) '.txt']));
end

disp('Done!');