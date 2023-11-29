%{
save the fitted OP for running MC simulaiton
generate random OP change according to the error PDF.
Also add more wavelength (outside the fitted wl)

Run 'S8_2_output_fitted_OP_random_PDF.m' before this script

Benjamin Kao
Last update: 2021/04/14
%}

clc;clear;close all;

global lambda net param_range 

%% param
to_output_wl=([650:1064])';

input_dir='fitted_result';
orig_output_dir='saved_OPs_random'; % the original output dir of 'S8_2_output_fitted_OP_random_PDF.m'
output_dir='saved_OPs_random_moreWL';
do_plot=1;
num_random_OP=30;

to_error_OP=[1 2 3 4 7 8]; % the OPs needs to add error, 1 = mua1, 2 = mus1, 3 = mua2 ......

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
all_sbj_OP_arr=[];
for i=1:length(target_name_arr)
    % load the previous saved subject OP information
    sbj_op_info=load(fullfile(input_dir,orig_output_dir,[target_name_arr{i} '_OP_info.mat']));
    
    % copy the previous saved subject OP information
    copyfile(fullfile(input_dir,orig_output_dir,[target_name_arr{i} '_OP_info.mat']),fullfile(input_dir,output_dir,[target_name_arr{i} '_OP_info.mat']));
    
    sbj_OP_arr=[];
    fitted_OP_index=1;
    
    choosed_rank=load(fullfile(input_dir,[target_name_arr{i} '_choosed_rank.txt']));
    
    % load the fitted tissue param
    sbj_fitting_result=load(fullfile(fitting_dir{target_fitting_dir_index(i)},'arrangement',[target_name_arr{i} '_fitting_result_sort.txt']));
    
    % load the OP error CDF IS
    OP_error_IS=load(fullfile(input_dir,[target_name_arr{i} '_error_IS_arr.mat']));
    OP_error_IS=OP_error_IS.error_IS_arr;
    
    % find the chosen fitted tissue param
    sbj_fitting_result=sbj_fitting_result(:,2:1+num_fitted_param); % delete the init index and the error
    choosed_fitting_result=sbj_fitting_result(choosed_rank,:);
    
    for j=1:length(choosed_rank)
        % make the random error
        random_OP_change_arr=sbj_op_info.OP_change_arr(2:end,:);
        
        % calculate the fitted OP
        fitted_mu=fun_param_to_mu(choosed_fitting_result(j,:),0);
        fitted_OP_arr=[to_output_wl fitted_mu];
        if length(all_sbj_OP_arr)==0
            all_sbj_OP_arr(:,:,1)=fitted_OP_arr;
        else
            all_sbj_OP_arr(:,:,end+1)=fitted_OP_arr;
        end
        sbj_OP_arr(:,:,fitted_OP_index)=fitted_OP_arr;
        save(fullfile(input_dir,output_dir,[target_name_arr{i} '_OP_' num2str(fitted_OP_index) '.txt' ]),'fitted_OP_arr','-ascii','-tabs');
        fitted_OP_index=fitted_OP_index+1;
        
        % add error to one OP
        for l=1:num_random_OP
            errored_OP_arr=fitted_OP_arr;
            errored_OP_arr(:,2:end)=errored_OP_arr(:,2:end).*random_OP_change_arr(l,:);
            sbj_OP_arr(:,:,fitted_OP_index)=errored_OP_arr;
            save(fullfile(input_dir,output_dir,[target_name_arr{i} '_OP_' num2str(fitted_OP_index) '.txt' ]),'errored_OP_arr','-ascii','-tabs');
            fitted_OP_index=fitted_OP_index+1;
        end
    end
    
    % plot the OPs
    if do_plot
        figure('Units','pixels','position',[0 0 1920 1080]);
        ti=tiledlayout('flow','TileSpacing','compact','Padding','none');
        OP_name_arr={'\mu_{a,scalp}','\mu_{s,scalp}','\mu_{a,skull}','\mu_{s,skull}','\mu_{a,CSF}','\mu_{s,CSF}','\mu_{a,GM}','\mu_{s,GM}',};
        for L=[1 2 3 4 7 8]
            nexttile();
            plot(sbj_OP_arr(:,1,1),squeeze(sbj_OP_arr(:,L+1,:)));
            title(OP_name_arr{L});
            set(gca,'fontsize',12, 'FontName', 'Times New Roman');
        end
        title(ti,strrep(target_name_arr{i},'_','\_'),'FontName', 'Times New Roman');
        print(fullfile(input_dir,output_dir,[target_name_arr{i} '_OP.png']),'-dpng','-r200');
        close all;
    end
end

%% calculate the average OP array
mean_OP=mean(all_sbj_OP_arr,3);

for i=1:length(target_name_arr)
    % save the mean OP
    save(fullfile(input_dir,output_dir,[target_name_arr{i} '_OP_' num2str(num_random_OP+2) '.txt' ]),'mean_OP','-ascii','-tabs');
end

disp('Done!');