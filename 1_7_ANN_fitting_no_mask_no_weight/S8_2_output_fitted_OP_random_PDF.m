%{
save the fitted OP for running MC simulaiton
generate random OP change according to the error PDF.

Benjamin Kao
Last update: 2021/03/18
%}

clc;clear;close all;

global lambda net param_range 

%% param
input_dir='fitted_result';
output_dir='saved_OPs_random';
do_plot=1;
num_random_OP=30; % number of random nosie to added to the OP

to_error_OP=[1 2 3 4 7 8]; % the OPs needs to add error, 1 = mua1, 2 = mus1, 3 = mua2 ......

target_name_arr={'tc_mean','ww_mean','wh_mean','yf_mean','kb_mean'}; % the name of the target spectrum

%% init
mkdir(fullfile(input_dir,output_dir));

%% main
all_sbj_OP_arr=[];
for i=1:length(target_name_arr)
    sbj_OP_arr=[];
    OP_change_arr=[]; % store the change (error) of each set of OP
    OP_answer_index_arr=[]; % store which answer is this OP based on
    fitted_OP_index=1;
    
    choosed_rank=load(fullfile(input_dir,[target_name_arr{i} '_choosed_rank.txt']));
    
    % load the OP error CDF IS
    OP_error_IS=load(fullfile(input_dir,[target_name_arr{i} '_error_IS_arr.mat']));
    OP_error_IS=OP_error_IS.error_IS_arr;
    
    for j=1:length(choosed_rank)
        % make the random error
        random_OP_change_arr=[];
        for opi=1:6
            random_OP_change_arr(:,opi)=1+interp1(linspace(0,1,length(OP_error_IS{opi})),OP_error_IS{opi},rand(num_random_OP,1),'pchip');
        end
        temp_arr=ones(num_random_OP,10);
        temp_arr(:,to_error_OP)=random_OP_change_arr;
        random_OP_change_arr=temp_arr;
        
        % load the fitted OP
        fitted_OP_arr=load(fullfile(input_dir,[target_name_arr{i} '_fitted_OP_' num2str(j) '.txt']));
        if length(all_sbj_OP_arr)==0
            all_sbj_OP_arr(:,:,1)=fitted_OP_arr;
        else
            all_sbj_OP_arr(:,:,end+1)=fitted_OP_arr;
        end
        sbj_OP_arr(:,:,fitted_OP_index)=fitted_OP_arr;
        save(fullfile(input_dir,output_dir,[target_name_arr{i} '_OP_' num2str(fitted_OP_index) '.txt' ]),'fitted_OP_arr','-ascii','-tabs');
        OP_change_arr(fitted_OP_index,:)=ones(1,10);
        OP_answer_index_arr(fitted_OP_index,1)=j;
        fitted_OP_index=fitted_OP_index+1;
        
        % add error to one OP
        for l=1:num_random_OP
            temp_OP_change_arr=random_OP_change_arr(l,:);
            OP_change_arr(fitted_OP_index,:)=temp_OP_change_arr;
            OP_answer_index_arr(fitted_OP_index,1)=j;
            errored_OP_arr=fitted_OP_arr;
            errored_OP_arr(:,2:end)=errored_OP_arr(:,2:end).*temp_OP_change_arr;
            sbj_OP_arr(:,:,fitted_OP_index)=errored_OP_arr;
            save(fullfile(input_dir,output_dir,[target_name_arr{i} '_OP_' num2str(fitted_OP_index) '.txt' ]),'errored_OP_arr','-ascii','-tabs');
            fitted_OP_index=fitted_OP_index+1;
        end
    end
    
    % save other information
    sbj_answer_number=length(choosed_rank);
    number_OP_set=size(sbj_OP_arr,3);
    save(fullfile(input_dir,output_dir,[target_name_arr{i} '_OP_info.mat']),'sbj_answer_number','number_OP_set','OP_change_arr','OP_answer_index_arr');
    
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
    % load the information for each subject
    load(fullfile(input_dir,output_dir,[target_name_arr{i} '_OP_info.mat']));
    number_OP_set=number_OP_set+1;
    OP_change_arr(end+1,:)=1;
    OP_answer_index_arr(end+1,:)=0;
    save(fullfile(input_dir,output_dir,[target_name_arr{i} '_OP_info.mat']),'sbj_answer_number','number_OP_set','OP_change_arr','OP_answer_index_arr');
    
    % save the mean OP
    save(fullfile(input_dir,output_dir,[target_name_arr{i} '_OP_' num2str(number_OP_set) '.txt' ]),'mean_OP','-ascii','-tabs');
end

disp('Done!');