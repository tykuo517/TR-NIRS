%{
Choose the fitted OP and output the result

Benjamin Kao
Last update: 2021/03/30
%}

clc;clear;close all; clearvars -global;

global lambda net param_range 

%% param
error_range_threshold=0.02; % The error difference smaller than this value will be consider multiple solution
output_dir='fitted_result'; % the folder to output the fitted OP
to_output_wl=(700:890)'; % The wavelength to output the fitted OP

fitting_dir={'fitting_SDS1234_3','fitting_SDS2345_3','fitting_SDS123456_3','fitting_SDS12345_3','fitting_SDS346_3'}; % The fitted folder for different SDS combination
fitted_OP_error=[13.2 5.7 26.3 10.7 23.1 56.2 % the error of fitted OP for each SDS combinations
                   15.4	8.1 26.4 14   24.5 51.7
                   12.4	5.4 24.7 10   19.8 40
                   12.5 5.4 26.5 10.5 22.9 44.9
                   15.7 7.2 30.4 14.0 22.4 55.8];
compensate_SDS_index_arr={[5 6],[1 6],[],[6],[1 2 5]}; % the SDSs didn't used to fitting

sbj_fitting_dir_arr=[3 5 1 4 1]; % which fitting folder does this subject use
target_name_arr={'tc_mean','ww_mean','wh_mean','yf_mean','kb_mean'}; % the name of the target spectrum
model_name_arr={'ZJ','WW2','WH2','YF','KB'}; % the name of ANN model corresponding to each target spec

toOutput_subject_index=1; % the index of the subject to output
% fix_toOutput_rank=[2 8]; % skip the chossing process, output these ranks
fix_toOutput_rank=[]; % skip the chossing process, output these ranks

num_SDS=6; % how many SDS are in the target spectrum
SDS_sim_correspond=[1 2 3 4 5 6]; % the SDS index in the simulated spectrum corresponding to each SDS in the target spectrum, SDS_sim_correspond(x)=y means 'measure x = sim y'

num_fitted_param=13; % the number of the fitted parameters

model_dir='model_arrange'; % the folder containing the arranged ANN file
target_dir='input_target_2'; % the folder containing the target spectrum

do_plot_anyPlot=1; % =1 to plot any plot
plot_shaded=0; % =1 to plot the fitted OP error

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
legend_height=50; % pixel, the height of legend box

plot_n_col=3; % the column number of subplot
plot_n_row=2; % the row number of subplot

%% init
fitting_index=sbj_fitting_dir_arr(toOutput_subject_index); % the fitting SDS sets (dir) used to fitting

if exist(output_dir,'dir')==0
    mkdir(output_dir);
end

fitting_wl=to_output_wl;
lambda=fitting_wl;
fun_init_param_to_mu_spec(); % load the epsilon for fitting wl

%% main
sbj=toOutput_subject_index;

fprintf('Process %s\n',target_name_arr{sbj});

% load ANN model
ANN_model=load(fullfile(model_dir,[model_name_arr{sbj} '_model.mat'])); % net, param_range
net=ANN_model.net; param_range=ANN_model.param_range;

% load the fitting result
fitting_result=load(fullfile(fitting_dir{fitting_index},'arrangement',[target_name_arr{sbj} '_fitting_result_sort.txt']));

% find the result ranks to output, also plot the find process
toOutput_rank=[]; % output these ranks

% plot
fig=figure('Units','pixels','position',[0 0 1000 500]);
set(fig,'visible','off');
ti=tiledlayout('flow','TileSpacing','compact','Padding','none');
nexttile(1);
hold on;
plot(1:size(fitting_result,1),fitting_result(:,end),'o','LineWidth',2);
xticks(1:size(fitting_result,1));
shadedErrorBar(1:size(fitting_result,1),fitting_result(1,end)*ones(1,size(fitting_result,1)),[1;0]*ones(1,size(fitting_result,1))*error_range_threshold,'lineProps',{'-','LineWidth',lineWidth},'patchSaturation',0.33);
xlabel('Rank');
ylabel('Fitting SDS error');
title('Compare fitting SDS error');
set(gca,'fontsize',12, 'FontName', 'Times New Roman');
grid on;

fprintf('Rank 1 error = %.2f%%\nRank 2 error = %.2f%%\n',fitting_result(1,end)*100,fitting_result(2,end)*100);
if fitting_result(1,end)<=fitting_result(2,end)-error_range_threshold
    fprintf('Rank 1 is smaller than threshold\n');
    toOutput_rank=1;
else
    in_threshold_index=find(fitting_result(:,end)<=fitting_result(1,end)+error_range_threshold);
    fprintf('Rank '); fprintf('%d ',in_threshold_index); fprintf('are in threshold.\n');
    if length(compensate_SDS_index_arr{fitting_index})==0
        toOutput_rank=in_threshold_index;
        fprintf('No other SDS, select '); fprintf('%d ',in_threshold_index); fprintf('as result\n');
    else
        compensate_SDS_error=fitting_result(in_threshold_index,1+num_fitted_param+compensate_SDS_index_arr{fitting_index});
        compensate_SDS_error_RMSE=sqrt(mean(compensate_SDS_error.^2,2));
        [compensate_SDS_error_RMSE_sort,compensate_order]=sort(compensate_SDS_error_RMSE,'ascend');
        fprintf('The error of other SDSs are '); fprintf('%.2f%% ',compensate_SDS_error_RMSE_sort*100); fprintf('\n');
        
        % plot
        nexttile(2);
        hold on;
        plot(in_threshold_index,compensate_SDS_error_RMSE,'o','LineWidth',2);
        xticks(in_threshold_index);
        shadedErrorBar(in_threshold_index,compensate_SDS_error_RMSE_sort(1)*ones(1,length(in_threshold_index)),[1;0]*ones(1,length(in_threshold_index))*error_range_threshold,'lineProps',{'-','LineWidth',lineWidth},'patchSaturation',0.33);
        xlabel('Rank');
        ylabel('Other SDS error');
        title('Compare other SDS error');
        set(gca,'fontsize',12, 'FontName', 'Times New Roman');
        grid on;
        
        if compensate_SDS_error_RMSE_sort(1)<=compensate_SDS_error_RMSE_sort(2)-error_range_threshold
            toOutput_rank=compensate_order(1);
            fprintf('Rank %d has other SDS error smaller than threshold\n',toOutput_rank);
        else
            in_threshold_index2=find(compensate_SDS_error_RMSE<=compensate_SDS_error_RMSE_sort(1)+error_range_threshold);
            toOutput_rank=sort(in_threshold_index2,'ascend');
            fprintf('Rank '); fprintf('%d ',toOutput_rank); fprintf('are in other SDS error threshold.\nSelect them as result\n');
        end
    end
end

print(fullfile(output_dir,[target_name_arr{sbj} '_choose_process.png']),'-dpng','-r200');
close all;

if length(fix_toOutput_rank)>=1
    toOutput_rank=fix_toOutput_rank;
end


% orig_init_index=fitting_result(:,1); % get the fitting index number
fitting_result=fitting_result(toOutput_rank,2:end); % delete the fitting index number

fitted_OP_arr=[];
fitted_spec_arr=[];
for i=1:size(fitting_result,1)
    temp_op=fun_param_to_mu(fitting_result(i,:),0);
    fitted_OP_arr(:,:,i)=temp_op;
    fitted_spec_arr(:,:,i)=fun_ANN_forward(temp_op);
end

%% save
for i=1:size(fitting_result,1)
    to_save=[lambda fitted_OP_arr(:,:,i)];
    save(fullfile(output_dir,[target_name_arr{sbj} '_fitted_OP_' num2str(i) '.txt']),'to_save','-ascii','-tabs');
    to_save=[lambda fitted_spec_arr(:,:,i)];
    save(fullfile(output_dir,[target_name_arr{sbj} '_fitted_spec_' num2str(i) '.txt']),'to_save','-ascii','-tabs');
end

save(fullfile(output_dir,[target_name_arr{sbj} '_choosed_rank.txt']),'toOutput_rank','-ascii','-tabs');
OP_CV=fitted_OP_error(fitting_index,:);
save(fullfile(output_dir,[target_name_arr{sbj} '_OP_CV.txt']),'OP_CV','-ascii','-tabs');
save(fullfile(output_dir,[target_name_arr{sbj} '_fitting_index.txt']),'fitting_index','-ascii','-tabs');
save(fullfile(output_dir,[target_name_arr{sbj} '_fitted_OP_info.mat']),'lambda','fitting_index','fitted_OP_arr','toOutput_rank','OP_CV');

%% plot the OPs
if do_plot_anyPlot
    fprintf('Plot fitted OP\n');
    legend_arr={};
    for j=1:length(toOutput_rank)
        legend_arr{end+1}=['rank ' num2str(toOutput_rank(j))]; % OP error
    end

    fig=figure('Units','pixels','position',[0 0 (left_spacing+subplot_width)*plot_n_col+right_spacing (upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing+lower_spacing]);
    set(fig,'visible','off');
    op_name={'\mu_{a','\mu_{s'};
    layer_name_arr={'scalp','skull','CSF','GM'};
    if length(toOutput_rank)>1
        colormap_arr=jet(length(toOutput_rank));
    else
        colormap_arr=jet(2);
        colormap_arr=colormap_arr(1,:);
    end
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
            for j=1:length(toOutput_rank)

                if plot_shaded
                    shadedErrorBar(fitting_wl,fitted_OP_arr(:,2*(L-1)+opi,j),fitted_OP_arr(:,2*(L-1)+opi,j)*fitted_OP_error(fitting_index,OP_index)*0.01,'lineProps',{'-','Color',colormap_arr(j,:),'LineWidth',lineWidth},'patchSaturation',0.33);
                else
                    plot(fitting_wl,fitted_OP_arr(:,2*(L-1)+opi,j),'-','Color',colormap_arr(j,:),'LineWidth',lineWidth);
                end
            end

            % scale to ub and lb
            ylim([param_range(2,(opi-1)*4+L) param_range(1,(opi-1)*4+L)]);
            xlabel('wavelength(nm)');
            ylabel([op_name{opi} ',' layer_name_arr{L} '}(1/cm)']);

%             if length(toOutput_rank)>1
%                 yyaxis right;
%                 CV_arr=std(squeeze(fitted_OP_arr(:,2*(L-1)+opi,:)),[],2)./mean(squeeze(fitted_OP_arr(:,2*(L-1)+opi,:)),2);
%                 plot(fitting_wl,CV_arr,'--');
%                 ylabel('CV');
%             end

%                         title(['layer ' num2str(L) ' ' op_name{opi}]);
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

    print(fullfile(output_dir,[target_name_arr{sbj} '_choosed_OP.png']),'-dpng','-r200');
    close all;
end

%% plot the fitted spec together
target_spec=load(fullfile(target_dir,[target_name_arr{sbj} '.txt']));

if do_plot_anyPlot
    fprintf('Plot the fitted spec\n');

    fig=figure('Units','pixels','position',[0 0 (left_spacing+subplot_width)*plot_n_col+right_spacing (upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing+lower_spacing]);
    set(fig,'visible','off');
    colormap_arr=lines(length(toOutput_rank)+1);
    % old version plot
%     legend_arr={'target'};
%     for j=1:length(toOutput_rank)
%         legend_arr{end+1}=['rank ' num2str(toOutput_rank(j)) ', spec error= ' num2str(fitting_result(j,num_fitted_param+num_SDS+1)*100,'%.2f%%')]; % spec error
%     end

%     subplot_index=1;
%     for s=1:num_SDS
%         row_index=ceil(subplot_index/plot_n_col);
%         col_index=subplot_index-(row_index-1)*plot_n_col;
%         subplot(2,3,subplot_index);
%         subplot_index=subplot_index+1;
%         hold on;
%         plot(target_spec(:,1),target_spec(:,s+1),'LineWidth',lineWidth);
%         for j=1:length(toOutput_rank)
%             plot(lambda,fitted_spec_arr(:,s,j),'--','Color',colormap_arr(j,:),'LineWidth',lineWidth);
%         end
%         xlabel('wavelength(nm)');
%         ylabel('reflectance');
%         title(['SDS ' num2str(s)]);
%         set(gca,'fontsize',fontSize, 'FontName', 'Times New Roman');
%         set(gca,'Unit','pixels','Position',[left_spacing+(left_spacing+subplot_width)*(col_index-1) lower_spacing+legend_height+upper_spacing+(subplot_height+upper_spacing)*(plot_n_row-row_index) subplot_width subplot_height]);
%         grid on;
%         if subplot_index==7
%             lgd=legend(legend_arr,'Location','southoutside','fontsize',lgdFontSize);
%             lgd.NumColumns=lgdNumCol;
%             set(lgd,'Unit','pixels','position',[left_spacing lower_spacing plot_n_col*subplot_width+(plot_n_col-1)*left_spacing legend_height]);
%         end
%     end
% 
%     % add the title
%     axes;
%     axis off;
%     set(gca,'Unit','pixels','Position',[0 0 (left_spacing+subplot_width)*plot_n_col+right_spacing (upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing+lower_spacing]);
%     text(((left_spacing+subplot_width)*plot_n_col+right_spacing)/2,(upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing*2/3+lower_spacing,[strrep(target_name_arr{sbj},'_',' ')],'Unit','pixels','fontsize',fontSize, 'FontName', 'Times New Roman','HorizontalAlignment','center','VerticalAlignment','middle')

    % new version plot
    SDS_cm_arr=[0.8 1.5 2.12 3 3.35 4.5];
    ti=tiledlayout('flow','TileSpacing','compact','Padding','none');
    for s=1:num_SDS
        nexttile();
        hold on;
        plot(target_spec(:,1),target_spec(:,s+1),'Color',colormap_arr(1,:),'LineWidth',lineWidth);
        for j=1:length(toOutput_rank)
            plot(lambda,fitted_spec_arr(:,s,j),'--','Color',colormap_arr(j+1,:),'LineWidth',lineWidth);
        end
        xlabel('wavelength(nm)');
        ylabel('reflectance');
        title(['SDS = ' num2str(SDS_cm_arr(s)) ' cm']);
        set(gca,'fontsize',fontSize, 'FontName', 'Times New Roman');
        grid on;
        yylim=ylim(); yylim(1)=0; ylim(yylim);
        
        legend_arr={'target'};
        for j=1:length(toOutput_rank)
            legend_arr{end+1}=['rank ' num2str(toOutput_rank(j)) ', spec error= ' num2str(fitting_result(j,num_fitted_param+s)*100,'%.2f%%')]; % spec error
        end
        lgd=legend(legend_arr,'Location','southoutside','fontsize',lgdFontSize);
        lgd.NumColumns=lgdNumCol;
    end

    print(fullfile(output_dir,[target_name_arr{sbj} '_choosed_spec.png']),'-dpng','-r200');
    close all;
end

disp('Done!');