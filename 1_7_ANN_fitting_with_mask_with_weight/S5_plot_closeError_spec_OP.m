%{
plot the fitting result close enough to the best result
For OP, also add the fitting noise for each OP on the plot

Benjamin Kao
Last update: 2021/01/19
%}

clc;clear;close all; clearvars -global;

global lambda net param_range 

%% param
error_range_threshold=0.02; % plot the fitting in this distance to the best result
max_choose_number=5; % the max number of fitting

fitting_dir={'fitting_SDS1234_3','fitting_SDS2345_3','fitting_SDS123456_3','fitting_SDS12345_3','fitting_SDS346_3'}; % please move the fitting folders into this folder first.
fitted_OP_error=[13.2 5.7 26.3 10.7 23.1 56.2 % the error of fitted OP
                   15.4	8.1 26.4 14   24.5 51.7
                   12.4	5.4 24.7 10   19.8 40
                   12.5 5.4 26.5 10.5 22.9 44.9
                   15.7 7.2 30.4 14.0 22.4 55.8];
fitting_index=3; % the fitting SDS sets (dir) used to fitting

target_name_arr={'tc_mean','ww_mean','wh_mean','yf_mean','kb_mean'}; % the name of the target spectrum
model_name_arr={'ZJ','WW2','WH2','YF','KB'}; % the name of ANN model corresponding to each target spec
times_to_fitting=20; % number of fitting using different init value

num_SDS=6; % how many SDS are in the target spectrum
SDS_sim_correspond=[1 2 3 4 5 6]; % the SDS index in the simulated spectrum corresponding to each SDS in the target spectrum, SDS_sim_correspond(x)=y means 'measure x = sim y'

num_fitted_param=13; % the number of the fitted parameters

fitting_wl_file='fitting_wl.txt'; % the file in the epsilon folder containing the fitting wavelength

model_dir='model_arrange'; % the folder containing the arranged ANN file
target_dir='input_target_2'; % the folder containing the target spectrum

do_plot_anyPlot=1; % =1 to plot any plot

fontSize=16;
lineWidth=2;
lgdFontSize=14;
lgdNumCol=10;

subplot_height=250; % pixel, the height of subplot
subplot_width=450; % pixel, the width of subplot
left_spacing=170; % pixel, the space between subplot and the left things
right_spacing=90; % pixel, the space right of the last column of subplot
upper_spacing=90; % pixel, the space between subplot and the upper things
lower_spacing=30; % pixel, the space below the legend box
legend_height=100; % pixel, the height of legend box

plot_n_col=3; % the column number of subplot
plot_n_row=2; % the row number of subplot

%% init
fitting_wl=load(fullfile('epsilon',fitting_wl_file));
lambda=fitting_wl;
fun_init_param_to_mu_spec(); % load the epsilon for fitting wl

%% main

for sbj=1:length(target_name_arr)
    % load ANN model
    ANN_model=load(fullfile(model_dir,[model_name_arr{sbj} '_model.mat'])); % net, param_range
    net=ANN_model.net; param_range=ANN_model.param_range;

    fprintf('Process %s\n',target_name_arr{sbj});

    target_spec=load(fullfile(target_dir,[target_name_arr{sbj} '.txt']));
    
    temp_target_folder=dir(fullfile(fitting_dir{fitting_index},[target_name_arr{sbj} '*']));
    assert(length(temp_target_folder)==1,['Error! more than 1 fitting folder for ' target_name_arr{sbj} '!']);
    temp_target_folder=fullfile(fitting_dir{fitting_index},temp_target_folder.name);
    
    fitting_result=load(fullfile(fitting_dir{fitting_index},'arrangement',[target_name_arr{sbj} '_fitting_result_sort.txt']));
    orig_init_index=fitting_result(:,1); % get the fitting index number
    fitting_result=fitting_result(:,2:end); % delete the fitting index number

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
        legend_arr={};
        for j=1:length(in_range_index)
            legend_arr{end+1}=['rank ' num2str(j)]; % OP error
        end

        fig=figure('Units','pixels','position',[0 0 (left_spacing+subplot_width)*plot_n_col+right_spacing (upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing+lower_spacing]);
        set(fig,'visible','off');
        op_name={'\mu_{a','\mu_{s'};
        layer_name_arr={'scalp','skull','CSF','GM'};
        colormap_arr=jet(length(in_range_index));
        subplot_index=1;
        for L=[1 2 4]
            for opi=1:2
                OP_index=(L-1)*2+opi;
                if L==4
                    OP_index=OP_index-2;
                end
                row_index=opi;
                col_index=L;
                if L==4
                    col_index=3;
                end
                axes;
                subplot_index=subplot_index+1;
                hold on;
                for j=1:length(in_range_index)
%                     plot(fitting_wl,fitted_OP_arr(:,2*(L-1)+opi,j),'--','Color',colormap_arr(j,:),'LineWidth',lineWidth);
                    shadedErrorBar(fitting_wl,fitted_OP_arr(:,2*(L-1)+opi,j),fitted_OP_arr(:,2*(L-1)+opi,j)*fitted_OP_error(fitting_index,OP_index)*0.01,'lineProps',{'-','Color',colormap_arr(j,:),'LineWidth',lineWidth},'patchSaturation',0.33);
                end

                % scale to ub and lb
                ylim([param_range(2,(opi-1)*4+L) param_range(1,(opi-1)*4+L)]);
                xlabel('wavelength(nm)');
                ylabel([op_name{opi} ',' layer_name_arr{L} '}(1/cm)']);

                yyaxis right;
                CV_arr=std(squeeze(fitted_OP_arr(:,2*(L-1)+opi,:)),[],2)./mean(squeeze(fitted_OP_arr(:,2*(L-1)+opi,:)),2);
                plot(fitting_wl,CV_arr,'--');
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
        text(((left_spacing+subplot_width)*plot_n_col+right_spacing)/2,(upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing*2/3+lower_spacing,[strrep(target_name_arr{sbj},'_',' ')],'Unit','pixels','fontsize',fontSize, 'FontName', 'Times New Roman','HorizontalAlignment','center','VerticalAlignment','middle');

        print(fullfile(fitting_dir{fitting_index},'arrangement',[target_name_arr{sbj} '_ranged_mu.png']),'-dpng','-r200');
        close all;

        %% plot the fitted spec together
        fprintf('Plot the fitted spec\n');
        legend_arr={'target'};
        for j=1:length(in_range_index)
            legend_arr{end+1}=['rank ' num2str(j) ', spec error= ' num2str(fitting_result(j,num_fitted_param+num_SDS+1)*100,'%.2f%%')]; % spec error
        end

        fig=figure('Units','pixels','position',[0 0 (left_spacing+subplot_width)*plot_n_col+right_spacing (upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing+lower_spacing]);
        set(fig,'visible','off');
        colormap_arr=jet(length(in_range_index));
        subplot_index=1;
        SDS_cm_arr=[0.8 1.5 2.12 3 3.35 4.5];
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
            title(['SDS = ' num2str(SDS_cm_arr(s)) ' cm']);
            set(gca,'fontsize',fontSize, 'FontName', 'Times New Roman');
            set(gca,'Unit','pixels','Position',[left_spacing+(left_spacing+subplot_width)*(col_index-1) lower_spacing+legend_height+upper_spacing+(subplot_height+upper_spacing)*(plot_n_row-row_index) subplot_width subplot_height]);
            grid on;
            if subplot_index==7
                lgd=legend(legend_arr,'Location','southoutside','fontsize',lgdFontSize);
                lgd.NumColumns=lgdNumCol;
                set(lgd,'Unit','pixels','position',[left_spacing lower_spacing plot_n_col*subplot_width+(plot_n_col-1)*left_spacing legend_height]);
            end
        end

%         % add the title
%         axes;
%         axis off;
%         set(gca,'Unit','pixels','Position',[0 0 (left_spacing+subplot_width)*plot_n_col+right_spacing (upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing+lower_spacing]);
%         text(((left_spacing+subplot_width)*plot_n_col+right_spacing)/2,(upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing*2/3+lower_spacing,[strrep(target_name_arr{sbj},'_',' ')],'Unit','pixels','fontsize',fontSize, 'FontName', 'Times New Roman','HorizontalAlignment','center','VerticalAlignment','middle')

        print(fullfile(fitting_dir{fitting_index},'arrangement',[target_name_arr{sbj} '_ranged_spec.png']),'-dpng','-r200');
        close all;
    end
end

disp('Done!');