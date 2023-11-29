%{
Calculate percentile of the best (smallest spec error) fitting result

Benjamin Kao
Last update: 2021/02/17
%}

clc;clear;close all; clearvars -global;

global lambda net param_range 

%% param
Add_error_mode=1; % if =0, use result without error; =1 to use results with error

input_dir='test_fitting_2023-10-31-21-57-22'; % please move the fitting folders into this folder first.
subject_name_arr={'KB'};
num_anser_to_generate=10; % number of target spec (true answer)
num_error_to_generate=15; % number of adding noise to the same, the first one will have no error
times_to_fitting=20; % number of fitting using different init value

fitting_dir='fitting_SDS123456_345_gate1-5';

% CW setting
num_SDS_cw=6; % how many SDS are in the target spectrum
SDS_choosed_cw=[1 2 3 4 5 6]; % the SDS chosen to fitted in the target spectrum
SDS_sim_correspond_cw=[1 2 3 4 5 6]; % the SDS index in the simulated spectrum corresponding to each SDS in the target spectrum, SDS_sim_correspond(x)=y means 'measure x = sim y'

% TR setting
num_SDS_tr=5;
SDS_choosed_tr=[]; % the SDS chosen to fitted in the target spectrum
SDS_sim_correspond_tr=[1 2 3 4 5]; % the SDS index in the simulated spectrum corresponding to each SDS in the target spectrum, SDS_sim_correspond(x)=y means 'measure x = sim y'

num_gate=10;
gate_choosed=[1 2 3 4 5];
gate_sim_correspond=[1 2 3 4 5 6 7 8 9 10];

% other setting
num_fitted_param=13; % the number of the fitted parameters

fitting_wl_file='fitting_wl.txt'; % the file in the epsilon folder containing the fitting wavelength

model_dir='model_arrange'; % the folder containing the arranged ANN file

do_plot_anyPlot=1; % =1 to plot any plot
do_plot_individual_fitting=1; % =1 to plot the result of each individual fitting

fontSize=14;
lineWidth=2;
lgdFontSize=14;
lgdNumCol=10;

subplot_height=300; % pixel, the height of subplot
subplot_width=450; % pixel, the width of subplot
left_spacing=120; % pixel, the space between subplot and the left things
right_spacing=50; % pixel, the space right of the last column of subplot
upper_spacing=150; % pixel, the space between subplot and the upper things
lower_spacing=30; % pixel, the space below the legend box
legend_height=0; % pixel, the height of legend box

plot_n_col=7; % the column number of subplot
plot_n_row=4; % the row number of subplot

%% init
if Add_error_mode==0
    to_process_fitting_index=1;
elseif Add_error_mode==1
    to_process_fitting_index=2:num_error_to_generate;
end

fitting_wl=load(fullfile('epsilon',fitting_wl_file));
lambda=fitting_wl;

mkdir(fullfile(input_dir,'arrangement',fitting_dir));

% load the param and OP answer
% param_answer_arr=load(fullfile(input_dir,'answers','param_answer_arr.txt'));
OP_answer_arr=[];
for i=1:num_anser_to_generate %1:num_anser_to_generate
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
error_percentile_arr=[]; % target * param * rank
OP_error=[];

for sbj=1:length(subject_name_arr)
    % load the fitting result
    for target_i=1:num_anser_to_generate %1:num_anser_to_generate
        fprintf('Process %s target %d\n',subject_name_arr{sbj},target_i);
        
        for error_i=to_process_fitting_index
            temp_target_folder=fullfile(input_dir,subject_name_arr{sbj},['target_' num2str(target_i) '_' num2str(error_i)],fitting_dir);
            
            temp_fitting_result=load(fullfile(temp_target_folder,'fitting_result_sort.txt'));
            temp_fitting_result=temp_fitting_result(:,2:end);
            
            temp_OP_error=temp_fitting_result(:,num_fitted_param+num_SDS_cw+num_SDS_tr+4:end);
            OP_error(end+1,:,:)=temp_OP_error;
            min_operr=min(temp_OP_error,[],1);
            max_operr=max(temp_OP_error,[],1);
            OP_error_arr=(temp_OP_error-min_operr)./(max_operr-min_operr);
            
            temp_percentile_arr=[];
            
            for i=1:times_to_fitting
                temp_percentile_arr(1,:,i)=OP_error_arr(i,:);
            end
            
            error_percentile_arr(end+1,:,:)=temp_percentile_arr;
        end
    end
end

%% save
mean_error_percentile=transpose(squeeze(mean(error_percentile_arr,1)));
std_error_percentile=transpose(squeeze(std(error_percentile_arr,[],1)));
if Add_error_mode==0
    save(fullfile(input_dir,'arrangement',fitting_dir,'opErrorPercentile_noError.mat'),'mean_error_percentile','std_error_percentile');
elseif Add_error_mode==1
    save(fullfile(input_dir,'arrangement',fitting_dir,'opErrorPercentile_withError.mat'),'mean_error_percentile','std_error_percentile');
end

%% plot
% plot percentile
figure('Position',[0 0 1920 1080]);
ti=tiledlayout('flow','TileSpacing','compact','Padding','none');
for i=1:times_to_fitting
    nexttile();
    hist(error_percentile_arr(:,4,i),20);
    title(['rank ' num2str(i)]);
    xlabel('normalized all OPs error');
    ylabel('count');
    xlim([0 1]);
    grid on;
    set(gca,'fontsize',fontSize, 'FontName', 'Times New Roman');
end

if Add_error_mode==0
    print(fullfile(input_dir,'arrangement',fitting_dir,'opErrorPercentile_noError.png'),'-dpng','-r200');
elseif Add_error_mode==1
    print(fullfile(input_dir,'arrangement',fitting_dir,'opErrorPercentile_withError.png'),'-dpng','-r200');
end

% plot OP error of best fitting result of each target 
title_arr={'\mu_{a,skull}','\mu_{s,skull}','\mu_{a,scalp}','\mu_{s,scalp}','\mu_{a,GM}','\mu_{s,GM}','all OPs'};
figure('Position',[0 0 1920 1080]);
ti=tiledlayout('flow','TileSpacing','compact','Padding','none');

OP_arr_for_yylim=OP_error(:,1,:);

for i=1:7
    nexttile();
    for j=1
        plot(100*OP_error(:,j,i),'o','LineWidth',2);
        hold on
    end
    title(title_arr(i));
    xlabel('number');
    ylabel('error(%)');
%     xlim([0 1]);
    ylim([min(100*OP_arr_for_yylim(:)) max(100*OP_arr_for_yylim(:))]);
%     xticks([1 2 3 4 5]);
    grid on;
    set(gca,'fontsize',fontSize, 'FontName', 'Times New Roman');
end

if Add_error_mode==0
    title(ti,'rank 1, no error');
    print(fullfile(input_dir,'arrangement',fitting_dir,'opError_noError.png'),'-dpng','-r200');
elseif Add_error_mode==1
    title(ti,'rank 1, with error');
    print(fullfile(input_dir,'arrangement',fitting_dir,'opError_withError.png'),'-dpng','-r200');
end

% close all;

disp('Done!');