%{
Save the literature OP for subjects to simulate

Benjamin Kao
Last update: 2021/03/18
%}

clc;clear;close all;

global lambda net param_range 

%% param
input_dir=fullfile('literature_OPs','OPs_to_sim_15');
output_dir=fullfile('fitted_result','saved_OPs_random');

target_name_arr={'tc_mean','ww_mean','wh_mean','yf_mean','kb_mean'}; % the name of the target spectrum

%% main
literature_wl=load(fullfile(input_dir,'sim_wl.txt'));
literature_op=load(fullfile(input_dir,'toSim_OP_1.txt'));
litOP=[literature_wl literature_op];

for i=1:length(target_name_arr)
    % load the information for each subject
    load(fullfile(output_dir,[target_name_arr{i} '_OP_info.mat']));
    number_OP_set=number_OP_set+1;
    OP_change_arr(end+1,:)=0;
    OP_answer_index_arr(end+1,:)=0;
    save(fullfile(output_dir,[target_name_arr{i} '_OP_info.mat']),'sbj_answer_number','number_OP_set','OP_change_arr','OP_answer_index_arr');
    
    % save the mean OP
    save(fullfile(output_dir,[target_name_arr{i} '_OP_' num2str(number_OP_set) '.txt' ]),'litOP','-ascii','-tabs');
end

disp('Done!');