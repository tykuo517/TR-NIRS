%{
Save the fitted OPs for wider WL, from the fitted tissue parameter to OPs

Benjamin Kao
Last update: 2021/02/25
%}

clc;clear;close all;

global lambda net param_range 

%% param
input_dir='fitted_result';
output_dir='saved_OPs_moreWL';
to_output_wl=(650:1000)';
do_plot=0;

num_fitted_param=13; % the number of the fitted parameters

fitting_dir={'fitting_SDS1234_3','fitting_SDS2345_3','fitting_SDS123456_3','fitting_SDS12345_3','fitting_SDS346_3'}; % please move the fitting folders into this folder first.
target_name_arr={'tc_mean','ww_mean','wh_mean','yf_mean','kb_mean'}; % the name of the target spectrum
target_fitting_dir_index=[3 5 1 4 1]; % the fitting dir for each target

%% init
if exist(fullfile(input_dir,output_dir),'dir')==0
    mkdir(fullfile(input_dir,output_dir));
end
lambda=to_output_wl;
fun_init_param_to_mu_spec(); % load the epsilon for fitting wl

%% main
for sbj=1:length(target_name_arr)
    % load which rank of the fitting result is chosen
    sbj_choosed_rank=load(fullfile(input_dir,[target_name_arr{sbj} '_choosed_rank.txt']));
    
    % load the fitted tissue param
    sbj_fitting_result=load(fullfile(fitting_dir{target_fitting_dir_index(sbj)},'arrangement',[target_name_arr{sbj} '_fitting_result_sort.txt']));
    
    % find the chosen fitted tissue param
    sbj_fitting_result=sbj_fitting_result(:,2:1+num_fitted_param); % delete the init index and the error
    choosed_fitting_result=sbj_fitting_result(sbj_choosed_rank,:);
    
    sbj_answer_number=length(sbj_choosed_rank);
    number_OP_set=length(sbj_choosed_rank);
    save(fullfile(input_dir,output_dir,[target_name_arr{sbj} '_OP_info.mat']),'sbj_answer_number','number_OP_set');
    
    % convert tissue param into OP
    for i=1:length(sbj_choosed_rank)
        fitted_mu=fun_param_to_mu(choosed_fitting_result(i,:),0);
        to_save=[to_output_wl fitted_mu];
        save(fullfile(input_dir,output_dir,[target_name_arr{sbj} '_OP_' num2str(i) '.txt']),'to_save','-ascii','-tabs');
        if do_plot
            % load the original output OPs
            orig_OP=load(fullfile(input_dir,[target_name_arr{sbj} '_fitted_OP_' num2str(i) '.txt']));
            
            figure('Units','pixels','position',[0 0 1920 1080]);
            ti=tiledlayout('flow','TileSpacing','compact','Padding','none');
            for j=1:10
                nexttile();
                plot(to_save(:,1),to_save(:,1+j),orig_OP(:,1),orig_OP(:,1+j));
            end
            title(ti,[strrep(target_name_arr{sbj},'_','\_') ' rank ' num2str(i)]);
            print(fullfile(input_dir,output_dir,[target_name_arr{sbj} '_OP_' num2str(i) '.png']),'-dpng','-r200');
            close all;
        end
    end
end

disp('Done!');