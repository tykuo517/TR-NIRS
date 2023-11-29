%{
Try to merge many fitting result, calculate the mean OP and forward the spec

Benjamin Kao
Last update: 2021/01/19
%}

clc;clear;close all; clearvars -global;

global lambda net param_range 

%% param
error_range_threshold=0.03; % plot the fitting in this distance to the best result
max_choose_number=5; % the max number of fitting

input_dir='test_fitting_2021-01-17-16-47-46'; % please move the fitting folders into this folder first.
subject_name_arr={'ZJ','WW2','WH2','YF','KB','YH','SJ','BT','SC'};
num_anser_to_generate=15; % number of target spec (true answer)
num_error_to_generate=15; % number of adding noise to the same, the first one will have no error
times_to_fitting=20; % number of fitting using different init value

fitting_dir='fitting_SDS1234';

num_SDS=6; % how many SDS are in the target spectrum
SDS_sim_correspond=[1 2 3 4 5 6]; % the SDS index in the simulated spectrum corresponding to each SDS in the target spectrum, SDS_sim_correspond(x)=y means 'measure x = sim y'

num_fitted_param=13; % the number of the fitted parameters

fitting_wl_file='fitting_wl.txt'; % the file in the epsilon folder containing the fitting wavelength

model_dir='model_arrange'; % the folder containing the arranged ANN file

do_plot_anyPlot=1; % =1 to plot any plot

fontSize=16;
lineWidth=2;
lgdFontSize=14;
lgdNumCol=10;

subplot_height=250; % pixel, the height of subplot
subplot_width=450; % pixel, the width of subplot
left_spacing=160; % pixel, the space between subplot and the left things
right_spacing=80; % pixel, the space right of the last column of subplot
upper_spacing=100; % pixel, the space between subplot and the upper things
lower_spacing=30; % pixel, the space below the legend box
legend_height=50; % pixel, the height of legend box

plot_n_col=4; % the column number of subplot
plot_n_row=3; % the row number of subplot

%% init
fitting_wl=load(fullfile('epsilon',fitting_wl_file));
lambda=fitting_wl;
fun_init_param_to_mu_spec(); % load the epsilon for fitting wl

mkdir(input_dir,'arrangement');

% load the param and OP answer
param_answer_arr=load(fullfile(input_dir,'answers','param_answer_arr.txt'));
OP_answer_arr=[];
for i=1:num_anser_to_generate
    temp_OP_arr=load(fullfile(input_dir,'answers',['OP_ans_' num2str(i) '.txt']));
    OP_answer_arr(:,:,i)=interp1(temp_OP_arr(:,1),temp_OP_arr(:,2:end),lambda);
end

%% main

% for sbj=1:length(subject_name_arr)
for sbj=[1]
    % make output folder and find fitting folder
    mkdir(fullfile(input_dir,'arrangement',subject_name_arr{sbj}));
    
    % load ANN model
    ANN_model=load(fullfile(model_dir,[subject_name_arr{sbj} '_model.mat'])); % net, param_range
    net=ANN_model.net; param_range=ANN_model.param_range;

    for target_i=1:num_anser_to_generate
        for error_i=1:num_error_to_generate
            fprintf('Process %s target %d, error %d\n',subject_name_arr{sbj},target_i,error_i);

            target_spec=load(fullfile(input_dir,subject_name_arr{sbj},['target_' num2str(target_i) '_' num2str(error_i)],'target_spec.txt'));
            target_spec_interp=interp1(target_spec(:,1),target_spec(:,2:end),lambda);
            
            temp_target_folder=fullfile(input_dir,subject_name_arr{sbj},['target_' num2str(target_i) '_' num2str(error_i)],fitting_dir);
            fitting_result=load(fullfile(temp_target_folder,'fitting_result_sort.txt'));
            orig_init_index=fitting_result(:,1);
            fitting_result=fitting_result(:,2:end);
            
            in_range_index=find(fitting_result(:,num_fitted_param+num_SDS+1)<=fitting_result(1,num_fitted_param+num_SDS+1)+error_range_threshold);
            assert(length(in_range_index)==max(in_range_index));
            if length(in_range_index)>max_choose_number
                in_range_index=in_range_index(1:max_choose_number);
            end
            
            fitted_spec_arr=[];
            fitted_OP_arr=[];
            for i=1:length(in_range_index)
                fitted_spec_arr(:,:,i)=load(fullfile(temp_target_folder,['fitting_' num2str(orig_init_index(i))],'fitted_spec.txt'));
                fitted_OP_arr(:,:,i)=load(fullfile(temp_target_folder,['fitting_' num2str(orig_init_index(i))],'fitted_mu.txt'));
            end
            
            mean_OP_arr=mean(fitted_OP_arr,3);
            mean_spec=fun_ANN_forward(mean_OP_arr);
            
%             figure;
%             for s=1:6
%                 subplot(2,3,s);
%                 plot(target_spec(:,1),target_spec(:,s+1),fitting_wl,mean_spec(:,s))
%             end
            

            if do_plot_anyPlot
                %% plot the mean OPs and spec
                fprintf('Plot fitted OP\n');

                fig=figure('Units','pixels','position',[0 0 (left_spacing+subplot_width)*plot_n_col+right_spacing (upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing+lower_spacing]);
                set(fig,'visible','off');
                op_name={'\mu_a','\mu_s'};
                subplot_index=1;
                for L=[1 2 4]
                    for opi=1:2
                        row_index=ceil(subplot_index/plot_n_col);
                        col_index=subplot_index-(row_index-1)*plot_n_col;
                        axes;
                        subplot_index=subplot_index+1;
                        hold on;
                        plot(fitting_wl,[OP_answer_arr(:,2*(L-1)+opi,target_i) mean_OP_arr(:,2*(L-1)+opi)],'LineWidth',lineWidth);
                        
                        % scale to ub and lb
                        ylim([param_range(2,(opi-1)*4+L) param_range(1,(opi-1)*4+L)]);
                        xlabel('wavelength(nm)');
                        ylabel([op_name{opi} '_' num2str(L) '(1/cm)']);
                        
                        yyaxis right;
                        error_arr=mean_OP_arr(:,2*(L-1)+opi)./OP_answer_arr(:,2*(L-1)+opi,target_i)-1;
                        plot(fitting_wl,error_arr*100);
                        ylabel('error(%)');
                        
                        title(['error = ' num2str(sqrt(mean(error_arr.^2))*100,'%.2f%%')]);
                        set(gca,'fontsize',fontSize, 'FontName', 'Times New Roman');
                        set(gca,'Unit','pixels','Position',[left_spacing+(left_spacing+subplot_width)*(col_index-1) lower_spacing+legend_height+upper_spacing+(subplot_height+upper_spacing)*(plot_n_row-row_index) subplot_width subplot_height]);
                        grid on;
                    end
                end
                
                fprintf('Plot fitted spec\n');
                
                for s=1:num_SDS
                    row_index=ceil(subplot_index/plot_n_col);
                    col_index=subplot_index-(row_index-1)*plot_n_col;
                    axes;
                    subplot_index=subplot_index+1;
                    hold on;
                    plot(target_spec(:,1),target_spec(:,s+1),'LineWidth',lineWidth);
                    plot(fitting_wl,mean_spec(:,s),'LineWidth',lineWidth);
                    xlabel('wavelength(nm)');
                    ylabel('reflectance');
                    
                    yyaxis right;
                    error_arr=mean_spec(:,s)./target_spec_interp(:,s)-1;
                    plot(fitting_wl,error_arr*100);
                    ylabel('error (%)');
                    
                    title(['SDS ' num2str(s) ', error= ' num2str(sqrt(mean(error_arr.^2))*100,'%.2f%%')]);
                    set(gca,'fontsize',fontSize, 'FontName', 'Times New Roman');
                    set(gca,'Unit','pixels','Position',[left_spacing+(left_spacing+subplot_width)*(col_index-1) lower_spacing+legend_height+upper_spacing+(subplot_height+upper_spacing)*(plot_n_row-row_index) subplot_width subplot_height]);
                    grid on;
                    
                    if s==num_SDS
                        lgd=legend({'answer','fitted'},'Location','southoutside','fontsize',lgdFontSize);
                        lgd.NumColumns=lgdNumCol;
                        set(lgd,'Unit','pixels','position',[left_spacing lower_spacing plot_n_col*subplot_width+(plot_n_col-1)*left_spacing legend_height]);
                    end
                end

                % add the title
                axes;
                axis off;
                set(gca,'Unit','pixels','Position',[0 0 (left_spacing+subplot_width)*plot_n_col+right_spacing (upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing+lower_spacing]);
                text(((left_spacing+subplot_width)*plot_n_col+right_spacing)/2,(upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing*2/3+lower_spacing,[subject_name_arr{sbj} ' target ' num2str(target_i) ' error ' num2str(error_i)],'Unit','pixels','fontsize',fontSize, 'FontName', 'Times New Roman','HorizontalAlignment','center','VerticalAlignment','middle');

%                 print(fullfile(temp_target_folder,'ranged_fitted_mu.png'),'-dpng','-r200');
                print(fullfile(input_dir,'arrangement',subject_name_arr{sbj},[subject_name_arr{sbj} '_' num2str(target_i) '_' num2str(error_i) '_tryMean.png']),'-dpng','-r200');
                close all;

            end
        end
    end
end

disp('Done!');