%{
Calculate the correlation coefficient of the OP error and the fitted spectrum error

Benjamin Kao
Last update: 2021/02/17
%}

clc;clear;close all; clearvars -global;

global lambda net param_range 

%% param
error_range_threshold=0.03; % plot the fitting in this distance to the best result
max_choose_number=5; % the max number of fitting

input_dir='test_fitting_2023-10-12-00-56-05'; % please move the fitting folders into this folder first.
subject_name_arr={'KB'};
num_anser_to_generate=15; % number of target spec (true answer)
num_error_to_generate=15; % number of adding noise to the same, the first one will have no error
times_to_fitting=20; % number of fitting using different init value

fitting_dir='fitting_SDS123456_345_gate1-5';

num_SDS=10; % how many SDS are in the target spectrum
SDS_sim_correspond=[1 2 3 4 5 6 7 8 9 10]; % the SDS index in the simulated spectrum corresponding to each SDS in the target spectrum, SDS_sim_correspond(x)=y means 'measure x = sim y'

num_fitted_param=5; % the number of the fitted parameters

fitting_wl_file='fitting_wl.txt'; % the file in the epsilon folder containing the fitting wavelength

model_dir='model_arrange'; % the folder containing the arranged ANN file

do_plot_anyPlot=0; % =1 to plot any plot
do_plot_individual_fitting=1; % =1 to plot the result of each individual fitting

fontSize=18;
lineWidth=2;
lgdFontSize=14;
lgdNumCol=10;

subplot_height=300; % pixel, the height of subplot
subplot_width=450; % pixel, the width of subplot
left_spacing=120; % pixel, the space between subplot and the left things
right_spacing=50; % pixel, the space right of the last column of subplot
upper_spacing=150; % pixel, the space between subplot and the upper things
lower_spacing=10; % pixel, the space below the legend box
legend_height=0; % pixel, the height of legend box

plot_n_col=7; % the column number of subplot
plot_n_row=4; % the row number of subplot

%% init
fitting_wl=load(fullfile('epsilon',fitting_wl_file));
lambda=fitting_wl;

mkdir(fullfile(input_dir,'arrangement',fitting_dir));

% load the param and OP answer
% param_answer_arr=load(fullfile(input_dir,'answers','param_answer_arr.txt'));
OP_answer_arr=[];
for i=1:num_anser_to_generate
    temp_OP_arr=load(fullfile(input_dir,'answers',['OP_ans_' num2str(i) '.txt']));
%     OP_answer_arr(:,:,i)=interp1(temp_OP_arr(:,1),temp_OP_arr(:,2:end),lambda);
    OP_answer_arr(:,:,i)=temp_OP_arr(:,2:end);
end

colormap_arr=jet(num_error_to_generate);

legend_arr={};
for error_i=1:num_error_to_generate
    legend_arr{error_i}=['error ' num2str(error_i)];
end

%% main
% for non adding error target
nonE_error_corr_coeff_arr=[]; % the corr coef of each error sets
nonE_ranged_error_corr_coeff_arr=[]; % the corr coef of the choosed error sets

% for added error target
error_corr_coeff_arr=[]; % the corr coef of each error sets
ranged_error_corr_coeff_arr=[]; % the corr coef of the choosed error sets

for sbj=1:length(subject_name_arr)
    % make output folder and find fitting folder
    
    % load ANN model
    ANN_model=load(fullfile(model_dir,[subject_name_arr{sbj} '_model.mat'])); % net, param_range
    net=ANN_model.net; param_range=ANN_model.param_range;

    % load the fitting result
    for target_i=1:num_anser_to_generate
        fprintf('Process %s target %d\n',subject_name_arr{sbj},target_i);
        
        fitted_result_arr=[];
        
        
        for error_i=1:num_error_to_generate

            temp_target_folder=fullfile(input_dir,subject_name_arr{sbj},['target_' num2str(target_i) '_' num2str(error_i)],fitting_dir);
            
            temp_fitting_result=load(fullfile(temp_target_folder,'fitting_result_sort.txt'));
            temp_fitting_result=temp_fitting_result(:,2:end);
            fitted_result_arr(:,:,error_i)=temp_fitting_result;
            
            % choose the closest error result
            in_range_index=find(temp_fitting_result(:,num_fitted_param+num_SDS+1)<=temp_fitting_result(1,num_fitted_param+num_SDS+1)+error_range_threshold);
            assert(length(in_range_index)==max(in_range_index));
            if length(in_range_index)>max_choose_number
                in_range_index=in_range_index(1:max_choose_number);
            end
            
            choosed_fitted_result=temp_fitting_result(in_range_index,:);
            
            temp_arr=zeros(1,4);
            temp_arr_2=zeros(1,4);
            
            for param_i=1:4
                temp=corrcoef(temp_fitting_result(:,num_fitted_param+num_SDS+1),temp_fitting_result(:,num_fitted_param+num_SDS+1+param_i));
                temp_arr(1,param_i)=temp(1,2);
                if length(in_range_index)>1
                    temp=corrcoef(choosed_fitted_result(:,num_fitted_param+num_SDS+1),choosed_fitted_result(:,num_fitted_param+num_SDS+1+param_i));
                    temp_arr_2(1,param_i)=temp(1,2);
                end
            end
            if error_i==1
                nonE_error_corr_coeff_arr(end+1,:)=temp_arr;
                if length(in_range_index)>1
                    nonE_ranged_error_corr_coeff_arr(end+1,:)=temp_arr_2;
                end
            else
                error_corr_coeff_arr(end+1,:)=temp_arr;
                if length(in_range_index)>1
                    ranged_error_corr_coeff_arr(end+1,:)=temp_arr_2;
                end
            end
        end
    end
end

%% plot
% disp('Plot the figure together.');
% fig=figure('Units','pixels','position',[0 0 (left_spacing+subplot_width)*plot_n_col+right_spacing (upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing+lower_spacing]);
% set(fig,'visible','off');
% op_name_arr={'\mu_{a,scalp}','\mu_{s,scalp}','\mu_{a,skull}','\mu_{s,skull}','\mu_{a,GM}','\mu_{s,GM}','all OPs'};
% % subplot_index=1;
% 
% % for not adding error, all fitting
% for i=1:4
%     row_index=1;
%     col_index=i;
%     axes;
%     hist(nonE_error_corr_coeff_arr(:,i),20);
%     grid on;
%     xlabel('r');
%     ylabel('count');
%     xlim([-1 1]);
%     title([op_name_arr{i} ' error']);
% %     title([op_name_arr{i} ', all fittings']);
%     set(gca,'fontsize',fontSize, 'FontName', 'Times New Roman');
%     set(gca,'Unit','pixels','Position',[left_spacing+(left_spacing+subplot_width)*(col_index-1) lower_spacing+legend_height+upper_spacing+(subplot_height+upper_spacing)*(plot_n_row-row_index) subplot_width subplot_height]);
% end
% 
% % for not adding error, in range fitting
% for i=1:4
%     row_index=2;
%     col_index=i;
%     axes;
%     hist(nonE_ranged_error_corr_coeff_arr(:,i),20);
%     grid on;
%     xlabel('r');
%     ylabel('count');
%     xlim([-1 1]);
% %     title([op_name_arr{i} ', better fittings']);
%     set(gca,'fontsize',fontSize, 'FontName', 'Times New Roman');
%     set(gca,'Unit','pixels','Position',[left_spacing+(left_spacing+subplot_width)*(col_index-1) lower_spacing+legend_height+upper_spacing+(subplot_height+upper_spacing)*(plot_n_row-row_index) subplot_width subplot_height]);
% end
% 
% % for adding error, all fitting
% for i=1:4
%     row_index=3;
%     col_index=i;
%     axes;
%     hist(error_corr_coeff_arr(:,i),20);
%     grid on;
%     xlabel('r');
%     ylabel('count');
%     xlim([-1 1]);
% %     title([op_name_arr{i} ', all fittings']);
%     set(gca,'fontsize',fontSize, 'FontName', 'Times New Roman');
%     set(gca,'Unit','pixels','Position',[left_spacing+(left_spacing+subplot_width)*(col_index-1) lower_spacing+legend_height+upper_spacing+(subplot_height+upper_spacing)*(plot_n_row-row_index) subplot_width subplot_height]);
% end
% 
% % for adding error, in range fitting
% for i=1:4
%     row_index=4;
%     col_index=i;
%     axes;
%     hist(ranged_error_corr_coeff_arr(:,i),20);
%     grid on;
%     xlabel('r');
%     ylabel('count');
%     xlim([-1 1]);
% %     title([op_name_arr{i} ', better fittings']);
%     set(gca,'fontsize',fontSize, 'FontName', 'Times New Roman');
%     set(gca,'Unit','pixels','Position',[left_spacing+(left_spacing+subplot_width)*(col_index-1) lower_spacing+legend_height+upper_spacing+(subplot_height+upper_spacing)*(plot_n_row-row_index) subplot_width subplot_height]);
% end
% 
% 
% % add the title
% axes;
% axis off;
% set(gca,'Unit','pixels','Position',[0 0 (left_spacing+subplot_width)*plot_n_col+right_spacing (upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing+lower_spacing]);
% text(((left_spacing+subplot_width)*plot_n_col+right_spacing)/2,(upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing*2/3+lower_spacing,'hist of corr coef','Unit','pixels','fontsize',fontSize, 'FontName', 'Times New Roman','HorizontalAlignment','center','VerticalAlignment','middle');
% 
% %                 print(fullfile(temp_target_folder,'ranged_fitted_mu.png'),'-dpng','-r200');
% print(fullfile(input_dir,'arrangement',fitting_dir,'corrcoef.png'),'-dpng','-r200');
% close all;

%% plot separate
disp('Plot the figure separate 1.');
op_name_arr={'\mu_{a,skull}','\mu_{a,GM}','\mu_{s,GM}','all OPs'};
plot_n_col=4; % the column number of subplot
plot_n_row=2; % the row number of subplot


% for not adding error, all fitting
fig=figure('Units','pixels','position',[0 0 (left_spacing+subplot_width)*plot_n_col+right_spacing (upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing+lower_spacing]);
% set(fig,'visible','off');
t1=tiledlayout("flow");
for i=1:4
    nexttile;
    hist(nonE_error_corr_coeff_arr(:,i),20);
    grid on;
    xlabel('r');
    ylabel('count');
    xlim([-1 1]);
    title(['r(' op_name_arr{i} ' error, spectral error)']);
    set(gca,'fontsize',fontSize, 'FontName', 'Times New Roman');
end

% add the title
title(t1,'hist of corr coef');
print(fullfile(input_dir,'arrangement',fitting_dir,'corrcoef_1.png'),'-dpng','-r200');
% close all;



%% for not adding error, in range fitting
disp('Plot the figure separate 2.');
fig=figure('Units','pixels','position',[0 0 (left_spacing+subplot_width)*plot_n_col+right_spacing (upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing+lower_spacing]);
t1=tiledlayout("flow");
% set(fig,'visible','off');
for i=1:4
    nexttile;
    hist(nonE_ranged_error_corr_coeff_arr(:,i),20);
    grid on;
    xlabel('r');
    ylabel('count');
    xlim([-1 1]);
    title(['r(' op_name_arr{i} ' error, spectral error)']);
    set(gca,'fontsize',fontSize, 'FontName', 'Times New Roman');
end

% add the title
title(t1,'hist of corr coef');
print(fullfile(input_dir,'arrangement',fitting_dir,'corrcoef_2.png'),'-dpng','-r200');

%% for adding error, all fitting
% disp('Plot the figure separate 3.');
% fig=figure('Units','pixels','position',[0 0 (left_spacing+subplot_width)*plot_n_col+right_spacing (upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing+lower_spacing]);
% set(fig,'visible','off');
% for i=1:7
%     row_index=ceil(i/plot_n_col);
%     col_index=i-(row_index-1)*plot_n_col;
%     axes;
%     hist(error_corr_coeff_arr(:,i),20);
%     grid on;
%     xlabel('r');
%     ylabel('count');
%     xlim([-1 1]);
%     title(['r(' op_name_arr{i} ' error, spectral error)']);
%     set(gca,'fontsize',fontSize, 'FontName', 'Times New Roman');
%     set(gca,'Unit','pixels','Position',[left_spacing+(left_spacing+subplot_width)*(col_index-1) lower_spacing+legend_height+upper_spacing+(subplot_height+upper_spacing)*(plot_n_row-row_index) subplot_width subplot_height]);
% end
% 
% % add the title
% axes;
% axis off;
% set(gca,'Unit','pixels','Position',[0 0 (left_spacing+subplot_width)*plot_n_col+right_spacing (upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing+lower_spacing]);
% text(((left_spacing+subplot_width)*plot_n_col+right_spacing)/2,(upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing*2/3+lower_spacing,'hist of corr coef','Unit','pixels','fontsize',fontSize, 'FontName', 'Times New Roman','HorizontalAlignment','center','VerticalAlignment','middle');
% print(fullfile(input_dir,'arrangement',fitting_dir,'corrcoef_3.png'),'-dpng','-r200');
% 
% %% for adding error, in range fitting
% disp('Plot the figure separate 4.');
% fig=figure('Units','pixels','position',[0 0 (left_spacing+subplot_width)*plot_n_col+right_spacing (upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing+lower_spacing]);
% set(fig,'visible','off');
% for i=1:7
%     row_index=ceil(i/plot_n_col);
%     col_index=i-(row_index-1)*plot_n_col;
%     axes;
%     hist(ranged_error_corr_coeff_arr(:,i),20);
%     grid on;
%     xlabel('r');
%     ylabel('count');
%     xlim([-1 1]);
%     title(['r(' op_name_arr{i} ' error, spectral error)']);
%     set(gca,'fontsize',fontSize, 'FontName', 'Times New Roman');
%     set(gca,'Unit','pixels','Position',[left_spacing+(left_spacing+subplot_width)*(col_index-1) lower_spacing+legend_height+upper_spacing+(subplot_height+upper_spacing)*(plot_n_row-row_index) subplot_width subplot_height]);
% end
% 
% % add the title
% axes;
% axis off;
% set(gca,'Unit','pixels','Position',[0 0 (left_spacing+subplot_width)*plot_n_col+right_spacing (upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing+lower_spacing]);
% text(((left_spacing+subplot_width)*plot_n_col+right_spacing)/2,(upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing*2/3+lower_spacing,'hist of corr coef','Unit','pixels','fontsize',fontSize, 'FontName', 'Times New Roman','HorizontalAlignment','center','VerticalAlignment','middle');
% print(fullfile(input_dir,'arrangement',fitting_dir,'corrcoef_4.png'),'-dpng','-r200');
% close all;

disp('Done!');