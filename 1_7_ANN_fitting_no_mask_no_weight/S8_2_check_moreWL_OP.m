%{
Check if the moreWL and original output OP are the same
run 'S8_2_output_fitted_OP_random_PDF.m' and 'S8_2_output_fitted_OP_random_PDF_moreWL.m' before theis script

Benjamin Kao
Last update: 2021/04/14
%}

clc;clear;close all;

%% param
input_dir_1=fullfile('fitted_result','saved_OPs_random');
input_dir_2=fullfile('fitted_result','saved_OPs_random_moreWL');
target_name_arr={'tc_mean','ww_mean','wh_mean','yf_mean','kb_mean'}; % the name of the target spectrum

%% main

figure('Units','pixels','position',[0 0 1920 1080]);
for sbj=1:length(target_name_arr)
    
%     for i=1:33
    for i=32
        OP_spec_1=load(fullfile(input_dir_1,[target_name_arr{sbj} '_OP_' num2str(i) '.txt']));
        OP_spec_2=load(fullfile(input_dir_2,[target_name_arr{sbj} '_OP_' num2str(i) '.txt']));
        
        for j=1:10
            subplot(2,5,j);
            plot(OP_spec_2(:,1),OP_spec_2(:,j+1),OP_spec_1(:,1),OP_spec_1(:,j+1));
        end
        drawnow;
        pause(0.3)
    end
end