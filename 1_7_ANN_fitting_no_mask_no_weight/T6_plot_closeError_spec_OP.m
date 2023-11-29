%{
plot the fitting result close enough to the best result

Benjamin Kao
Last update: 2021/03/30
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
SDS_sim_correspond=[1 2 3 4]; % the SDS index in the simulated spectrum corresponding to each SDS in the target spectrum, SDS_sim_correspond(x)=y means 'measure x = sim y'
SDS_length_arr=[0.8 1.5 2.12 3 3.35 4.5]; % cm

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
upper_spacing=90; % pixel, the space between subplot and the upper things
lower_spacing=30; % pixel, the space below the legend box
legend_height=100; % pixel, the height of legend box

plot_n_col=3; % the column number of subplot
plot_n_row=2; % the row number of subplot

%% init
fitting_wl=load(fullfile('epsilon',fitting_wl_file));
lambda=fitting_wl;
fun_init_param_to_mu_spec(); % load the epsilon for fitting wl

mkdir(fullfile(input_dir,'arrangement',fitting_dir));

% load the param and OP answer
param_answer_arr=load(fullfile(input_dir,'answers','param_answer_arr.txt'));
OP_answer_arr=[];
for i=1:num_anser_to_generate
    temp_OP_arr=load(fullfile(input_dir,'answers',['OP_ans_' num2str(i) '.txt']));
    OP_answer_arr(:,:,i)=interp1(temp_OP_arr(:,1),temp_OP_arr(:,2:end),lambda);
end

%% main

for sbj=1:length(subject_name_arr)
    % make output folder and find fitting folder
    mkdir(fullfile(input_dir,'arrangement',fitting_dir,subject_name_arr{sbj}));
    
    % load ANN model
    ANN_model=load(fullfile(model_dir,[subject_name_arr{sbj} '_model.mat'])); % net, param_range
    net=ANN_model.net; param_range=ANN_model.param_range;

    for target_i=1:num_anser_to_generate
        for error_i=1:num_error_to_generate
            fprintf('Process %s target %d, error %d\n',subject_name_arr{sbj},target_i,error_i);

            target_spec=load(fullfile(input_dir,subject_name_arr{sbj},['target_' num2str(target_i) '_' num2str(error_i)],'target_spec.txt'));
            
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

            if do_plot_anyPlot
                %% plot the OPs
                fprintf('Plot fitted OP\n');
                legend_arr={'answer'};
                for j=1:length(in_range_index)
                    legend_arr{end+1}=['rank ' num2str(j) ', OP error= ' num2str(fitting_result(j,end)*100,'%.2f%%')]; % OP error
                end

                fig=figure('Units','pixels','position',[0 0 (left_spacing+subplot_width)*plot_n_col+right_spacing (upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing+lower_spacing]);
                set(fig,'visible','off');
                op_name={'\mu_{a','\mu_{s'};
                layer_name_arr={'scalp','skull','CSF','GM'};
                colormap_arr=jet(length(in_range_index));
                subplot_index=1;
                for L=[1 2 4]
                    for opi=1:2
                        row_index=opi;
                        col_index=L;
                        if col_index==4
                            col_index=3;
                        end
                        axes();
                        subplot_index=subplot_index+1;
                        hold on;
                        plot(fitting_wl,OP_answer_arr(:,2*(L-1)+opi,target_i),'LineWidth',lineWidth);
                        for j=1:length(in_range_index)
                            plot(fitting_wl,fitted_OP_arr(:,2*(L-1)+opi,j),'--','Color',colormap_arr(j,:),'LineWidth',lineWidth);
                        end
                        
                        % scale to ub and lb
                        ylim([param_range(2,(opi-1)*4+L) param_range(1,(opi-1)*4+L)]);
                        xlabel('wavelength(nm)');
                        ylabel([op_name{opi} ',' layer_name_arr{L} '}(1/cm)']);
                        
                        yyaxis right;
                        CV_arr=std(squeeze(fitted_OP_arr(:,2*(L-1)+opi,:)),[],2)./mean(squeeze(fitted_OP_arr(:,2*(L-1)+opi,:)),2);
                        plot(fitting_wl,CV_arr);
                        ylabel('CV');
                        
                        set(gca,'fontsize',fontSize, 'FontName', 'Times New Roman');
                        set(gca,'Unit','pixels','Position',[left_spacing+(left_spacing+subplot_width)*(col_index-1) lower_spacing+legend_height+upper_spacing+(subplot_height+upper_spacing)*(plot_n_row-row_index) subplot_width subplot_height]);
                        grid on;
                        if subplot_index==7
                            lgd=legend(legend_arr,'Location','southoutside','fontsize',lgdFontSize);
                            lgd.NumColumns=lgdNumCol;
                            set(lgd,'Unit','pixels','position',[left_spacing lower_spacing plot_n_col*subplot_width+(plot_n_col-1)*left_spacing legend_height]);
                        end
                    end
                end

                % add the title
                axes;
                axis off;
                set(gca,'Unit','pixels','Position',[0 0 (left_spacing+subplot_width)*plot_n_col+right_spacing (upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing+lower_spacing]);
                text(((left_spacing+subplot_width)*plot_n_col+right_spacing)/2,(upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing*2/3+lower_spacing,[subject_name_arr{sbj} ' target ' num2str(target_i) ' error ' num2str(error_i)],'Unit','pixels','fontsize',fontSize, 'FontName', 'Times New Roman','HorizontalAlignment','center','VerticalAlignment','middle');

%                 print(fullfile(temp_target_folder,'ranged_fitted_mu.png'),'-dpng','-r200');
                print(fullfile(input_dir,'arrangement',fitting_dir,subject_name_arr{sbj},[subject_name_arr{sbj} '_' num2str(target_i) '_' num2str(error_i) '_ranged_mu.png']),'-dpng','-r200');
                close all;

                %% plot the fitted spec together
                fprintf('Plot the fitted spec\n');
                legend_arr={'target'};
                for j=1:length(in_range_index)
                    legend_arr{end+1}=['rank ' num2str(j) ', spec error= ' num2str(fitting_result(j,num_fitted_param+num_SDS+1)*100,'%.2f%%')]; % spec error
                end

                fig=figure('Units','pixels','position',[0 0 (left_spacing+subplot_width)*plot_n_col+right_spacing (upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing+lower_spacing]);
                set(fig,'visible','off');
                op_name={'mua','mus'};
                colormap_arr=jet(length(in_range_index));
                subplot_index=1;
                for s=1:num_SDS
                    row_index=ceil(subplot_index/plot_n_col);
                    col_index=subplot_index-(row_index-1)*plot_n_col;
                    subplot(2,3,subplot_index);
                    subplot_index=subplot_index+1;
                    hold on;
                    plot(target_spec(:,1),target_spec(:,s+1),'LineWidth',lineWidth);
                    for j=1:length(in_range_index)
                        plot(lambda,fitted_spec_arr(:,s,j),'--','Color',colormap_arr(j,:),'LineWidth',lineWidth);
                    end
                    xlabel('wavelength(nm)');
                    ylabel('reflectance');
                    title(['SDS ' num2str(SDS_length_arr(s)) ' cm']);
                    set(gca,'fontsize',fontSize, 'FontName', 'Times New Roman');
                    set(gca,'Unit','pixels','Position',[left_spacing+(left_spacing+subplot_width)*(col_index-1) lower_spacing+legend_height+upper_spacing+(subplot_height+upper_spacing)*(plot_n_row-row_index) subplot_width subplot_height]);
                    grid on;
                    if subplot_index==7
                        lgd=legend(legend_arr,'Location','southoutside','fontsize',lgdFontSize);
                        lgd.NumColumns=lgdNumCol;
                        set(lgd,'Unit','pixels','position',[left_spacing lower_spacing plot_n_col*subplot_width+(plot_n_col-1)*left_spacing legend_height]);
                    end
                end

                % add the title
%                 axes;
%                 axis off;
%                 set(gca,'Unit','pixels','Position',[0 0 (left_spacing+subplot_width)*plot_n_col+right_spacing (upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing+lower_spacing]);
%                 text(((left_spacing+subplot_width)*plot_n_col+right_spacing)/2,(upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing*2/3+lower_spacing,[subject_name_arr{sbj} ' target ' num2str(target_i) ' error ' num2str(error_i)],'Unit','pixels','fontsize',fontSize, 'FontName', 'Times New Roman','HorizontalAlignment','center','VerticalAlignment','middle')

%                 print(fullfile(temp_target_folder,'ranged_fitted_spec.png'),'-dpng','-r200');
                print(fullfile(input_dir,'arrangement',fitting_dir,subject_name_arr{sbj},[subject_name_arr{sbj} '_' num2str(target_i) '_' num2str(error_i) '_ranged_spec.png']),'-dpng','-r200');
                close all;
            end
            
            %% plot the fitted error and the init error, also the OP error and the fitted error
            if do_plot_anyPlot
                fprintf('Plot the fitted param\n');
                fig=figure('Units','pixels','position',[0 0 1600 1200]);
                set(fig,'visible','off');
                ti=tiledlayout('flow','TileSpacing','compact','Padding','none');
                % error for each param
                param_name_arr={'hc1','sto2 1','hc2','sto2 2','hc4','sto2 4','mel'};
                for param_i=7:num_fitted_param % only plot the mua param
                    nexttile();
                    plot(fitting_result(1:length(in_range_index),num_fitted_param+num_SDS+1),fitting_result(1:length(in_range_index),param_i),'o','LineWidth',lineWidth);
                    grid on;
                    xlabel('spectral error');
                    ylabel('percentage error');
                    yline(param_answer_arr(target_i,param_i),'b','LineWidth',lineWidth);
                    title(param_name_arr{param_i-6});
                    set(gca,'fontsize',fontSize, 'FontName', 'Times New Roman');
                end
                % OP error
                param_name_arr={'\mu_{a,scalp}','\mu_{s,scalp}','\mu_{a,skull}','\mu_{s,skull}','\mu_{a,GM}','\mu_{s,GM}','all OPs'};
                for param_i=1:7
                    nexttile();
                    plot(fitting_result(1:length(in_range_index),num_fitted_param+num_SDS+1),fitting_result(1:length(in_range_index),num_fitted_param+num_SDS+1+param_i),'o','LineWidth',lineWidth);
                    grid on;
                    xlabel('spectral error');
                    ylabel('RMSPE');
                    title(param_name_arr{param_i});
                    yylim=ylim(); yylim(1)=0; ylim(yylim);
                    set(gca,'fontsize',fontSize, 'FontName', 'Times New Roman');
                end
%                 print(fullfile(temp_target_folder,'ranged_fitted_param.png'),'-dpng','-r200');
                print(fullfile(input_dir,'arrangement',fitting_dir,subject_name_arr{sbj},[subject_name_arr{sbj} '_' num2str(target_i) '_' num2str(error_i) '_ranged_param.png']),'-dpng','-r200');
                close all;
            end
        end
    end
end

disp('Done!');