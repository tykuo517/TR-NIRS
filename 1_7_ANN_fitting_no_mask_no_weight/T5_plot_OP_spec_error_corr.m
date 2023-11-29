%{
plot the correlation of the OP error and the fitted spectrum error

Benjamin Kao
Last update: 2021/01/19
%}

clc;clear;close all; clearvars -global;

global lambda net param_range 

%% param
input_dir='test_fitting_2021-01-22-03-42-40'; % please move the fitting folders into this folder first.
subject_name_arr={'ZJ','WW2','WH2','YF','KB','YH','SJ','BT','SC'};
num_anser_to_generate=15; % number of target spec (true answer)
num_error_to_generate=15; % number of adding noise to the same, the first one will have no error
times_to_fitting=20; % number of fitting using different init value

fitting_dir='fitting_SDS346';

num_SDS=6; % how many SDS are in the target spectrum
SDS_sim_correspond=[1 2 3 4 5 6]; % the SDS index in the simulated spectrum corresponding to each SDS in the target spectrum, SDS_sim_correspond(x)=y means 'measure x = sim y'

num_fitted_param=13; % the number of the fitted parameters

fitting_wl_file='fitting_wl.txt'; % the file in the epsilon folder containing the fitting wavelength

model_dir='model_arrange'; % the folder containing the arranged ANN file

do_plot_anyPlot=0; % =1 to plot any plot
do_plot_individual_fitting=1; % =1 to plot the result of each individual fitting

fontSize=16;
lineWidth=2;
lgdFontSize=14;
lgdNumCol=10;

subplot_height=250; % pixel, the height of subplot
subplot_width=450; % pixel, the width of subplot
left_spacing=80; % pixel, the space between subplot and the left things
right_spacing=50; % pixel, the space right of the last column of subplot
upper_spacing=90; % pixel, the space between subplot and the upper things
lower_spacing=30; % pixel, the space below the legend box
legend_height=100; % pixel, the height of legend box

plot_n_col=4; % the column number of subplot
plot_n_row=2; % the row number of subplot

%% init
fitting_wl=load(fullfile('epsilon',fitting_wl_file));
lambda=fitting_wl;

mkdir(fullfile(input_dir,'arrangement',fitting_dir));

% load the param and OP answer
param_answer_arr=load(fullfile(input_dir,'answers','param_answer_arr.txt'));
OP_answer_arr=[];
for i=1:num_anser_to_generate
    temp_OP_arr=load(fullfile(input_dir,'answers',['OP_ans_' num2str(i) '.txt']));
    OP_answer_arr(:,:,i)=interp1(temp_OP_arr(:,1),temp_OP_arr(:,2:end),lambda);
end

colormap_arr=jet(num_error_to_generate);

legend_arr={};
for error_i=1:num_error_to_generate
    legend_arr{error_i}=['error ' num2str(error_i)];
end

%% main

for sbj=1:length(subject_name_arr)
    % make output folder and find fitting folder
    
    % load ANN model
    ANN_model=load(fullfile(model_dir,[subject_name_arr{sbj} '_model.mat'])); % net, param_range
    net=ANN_model.net; param_range=ANN_model.param_range;

    % load the fitting result
%     for target_i=1:num_anser_to_generate
    for target_i=1
        fprintf('Process %s target %d\n',subject_name_arr{sbj},target_i);
        
        fitted_result_arr=[];
        
        for error_i=1:num_error_to_generate

            temp_target_folder=fullfile(input_dir,subject_name_arr{sbj},['target_' num2str(target_i) '_' num2str(error_i)],fitting_dir);
            
            temp_fitting_result=load(fullfile(temp_target_folder,'fitting_result_sort.txt'));
            fitted_result_arr(:,:,error_i)=temp_fitting_result(:,2:end);
        end
        
        % plot 
        fig=figure('Units','pixels','position',[0 0 (left_spacing+subplot_width)*plot_n_col+right_spacing (upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing+lower_spacing]);
        set(fig,'visible','off');
        ti=tiledlayout('flow','TileSpacing','compact','Padding','none');
        param_name_arr={'mua 1','mus 1','mua 2','mus 2','mua 4','mus 4','all OP'};
        for i=1:7
            row_index=ceil(i/plot_n_col);
            col_index=i-(row_index-1)*plot_n_col;
            axes;
%             nexttile();
            hold on;
            for error_i=1:num_error_to_generate
                plot(fitted_result_arr(:,num_fitted_param+num_SDS+1,error_i),fitted_result_arr(:,num_fitted_param+num_SDS+1+i,error_i),':o','Color',colormap_arr(error_i,:));
            end
            xlabel('fitted error');
            ylabel('error');
            title(param_name_arr{i});
            set(gca,'fontsize',fontSize, 'FontName', 'Times New Roman');
            grid on;
            set(gca,'Unit','pixels','Position',[left_spacing+(left_spacing+subplot_width)*(col_index-1) lower_spacing+legend_height+upper_spacing+(subplot_height+upper_spacing)*(plot_n_row-row_index) subplot_width subplot_height]);
            if i==7
                lgd=legend(legend_arr,'Location','southoutside','fontsize',lgdFontSize);
                lgd.NumColumns=lgdNumCol;
                set(lgd,'Unit','pixels','position',[left_spacing lower_spacing plot_n_col*subplot_width+(plot_n_col-1)*left_spacing legend_height]);
            end
        end
        
        % add the title
        axes;
        axis off;
        set(gca,'Unit','pixels','Position',[0 0 (left_spacing+subplot_width)*plot_n_col+right_spacing (upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing+lower_spacing]);
        text(((left_spacing+subplot_width)*plot_n_col+right_spacing)/2,(upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing*2/3+lower_spacing,[subject_name_arr{sbj} ' target ' num2str(target_i)],'Unit','pixels','fontsize',fontSize, 'FontName', 'Times New Roman','HorizontalAlignment','center','VerticalAlignment','middle')
        
        print(fullfile(input_dir,'arrangement',fitting_dir,[subject_name_arr{sbj} '_target_' num2str(target_i) '.png']),'-dpng','-r200');
        close all;
    end
end

disp('Done!');