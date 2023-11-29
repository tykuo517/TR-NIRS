%{
save the fitted OP for running MC simulaiton, also calculate the mean OP and plot, compare them

Benjamin Kao
Last update: 2021/03/17
%}

clc;clear;close all;

global lambda net param_range 

%% param
input_dir='fitted_result';
output_dir='saved_OPs';
do_plot=1;

add_OP_error=1; % if =1, add the OP error to the fitted result
to_error_OP=[1 2 3 4 7 8]; % the OPs needs to add error, 1 = mua1, 2 = mus1, 3 = mua2 ......

target_name_arr={'tc_mean','ww_mean','wh_mean','yf_mean','kb_mean'}; % the name of the target spectrum

%% init
mkdir(fullfile(input_dir,output_dir));

%% main
all_sbj_OP_arr=[];
all_answer_index_arr=[];
all_answer_subject_arr=[];

for i=1:length(target_name_arr)
    OP_answer_index_arr=[]; % store which answer is this OP based on
    fitted_OP_index=1;
    
    % load the fitted OP
    choosed_rank=load(fullfile(input_dir,[target_name_arr{i} '_choosed_rank.txt']));
    for j=1:length(choosed_rank)
        fitted_OP_arr=load(fullfile(input_dir,[target_name_arr{i} '_fitted_OP_' num2str(j) '.txt']));
        if length(all_sbj_OP_arr)==0
            all_sbj_OP_arr(:,:,1)=fitted_OP_arr;
        else
            all_sbj_OP_arr(:,:,end+1)=fitted_OP_arr;
        end
        save(fullfile(input_dir,output_dir,[target_name_arr{i} '_OP_' num2str(fitted_OP_index) '.txt' ]),'fitted_OP_arr','-ascii','-tabs');
        OP_answer_index_arr(fitted_OP_index,1)=j;
        all_answer_subject_arr(end+1,1)=i;
        all_answer_index_arr(end+1,1)=j;
        fitted_OP_index=fitted_OP_index+1;
    end
    
    % save other information
    sbj_answer_number=length(choosed_rank);
    save(fullfile(input_dir,output_dir,[target_name_arr{i} '_OP_info.mat']),'sbj_answer_number','OP_answer_index_arr');
end

%% calculate the average OP array
mean_OP=mean(all_sbj_OP_arr,3);

legend_arr={'mean'};
for i=1:size(all_sbj_OP_arr,3)
    if all_answer_index_arr(i)==1
        legend_arr{end+1}=['subject ' num2str(all_answer_subject_arr(i))];
    else
        legend_arr{end+1}=['subject ' num2str(all_answer_subject_arr(i)) '_' num2str(all_answer_index_arr(i))];
    end
end

if do_plot
    figure('Units','pixels','position',[0 0 1920 1080]);
    ti=tiledlayout('flow','TileSpacing','compact','Padding','none');
    OP_name_arr={'\mu_{a,scalp}','\mu_{s,scalp}','\mu_{a,skull}','\mu_{s,skull}','\mu_{a,CSF}','\mu_{s,CSF}','\mu_{a,GM}','\mu_{s,GM}',};
    for L=[1 2 3 4 7 8]
        nexttile();
        plot(mean_OP(:,1),mean_OP(:,L+1),'LineWidth',2);
        hold on;
        for i=1:size(all_sbj_OP_arr,3)
            plot(all_sbj_OP_arr(:,1,i),all_sbj_OP_arr(:,L+1,i));
        end
        title(OP_name_arr{L});
        xlabel('wavelength(nm)');
        ylabel([OP_name_arr{L} '(1/cm)']);
        lgd=legend(legend_arr,'Location','southoutside');
        lgd.NumColumns = 10;
        set(gca,'fontsize',12, 'FontName', 'Times New Roman');
    end
    title(ti,'Mean OP','FontName', 'Times New Roman');
    print(fullfile(input_dir,output_dir,'mean_OP.png'),'-dpng','-r200');
    close all;
end

for i=1:length(target_name_arr)
    % load the information for each subject
    load(fullfile(input_dir,output_dir,[target_name_arr{i} '_OP_info.mat']));
    OP_answer_index_arr(end+1,:)=0;
    save(fullfile(input_dir,output_dir,[target_name_arr{i} '_OP_info.mat']),'sbj_answer_number','OP_answer_index_arr');
    
    % save the mean OP
    save(fullfile(input_dir,output_dir,[target_name_arr{i} '_OP_' num2str(length(OP_answer_index_arr)) '.txt' ]),'mean_OP','-ascii','-tabs');
end

disp('Done!');