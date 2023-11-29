%{
Arrange the fitting result

Benjamin Kao
Last update: 2021/01/25
%}

clc;clear;close all; clearvars -global;

global lambda net param_range 

%% param
input_dir='fitting_SDS1234_3'; % The fitting dir

target_name_arr={'tc_mean','ww_mean','wh_mean','yf_mean','kb_mean'}; % the name of the target spectrum
model_name_arr={'ZJ','WW2','WH2','YF','KB'}; % the name of ANN model corresponding to each target spec

% SDS_choosed=[1 2 3 4 5 6]; % the SDS chosen to fitted in the target spectrum
SDS_choosed=[1 2 3 4]; % the SDS chosen to fitted in the target spectrum

num_SDS=6; % how many SDS are in the target spectrum
SDS_dist_arr=[0.8 1.5 2.12 3 3.35 4.5 4.74]; % cm
SDS_sim_correspond=[1 2 3 4 5 6]; % the SDS index in the simulated spectrum corresponding to each SDS in the target spectrum, SDS_sim_correspond(x)=y means 'measure x = sim y'

num_fitted_param=13; % the number of the fitted parameters
times_to_fitting=20; % number of init value used to fitting

model_dir='model_arrange'; % the folder containing the arranged ANN file
target_dir='input_target_2'; % the folder containing the target spectrum

do_plot_anyPlot=1; % =1 to plot any plot
do_plot_individual_fitting=1; % =1 to plot the result of each individual fitting

fontSize=16;
lineWidth=2;
lgdFontSize=14;
lgdNumCol=10;

subplot_height=250; % pixel, the height of subplot
subplot_width=450; % pixel, the width of subplot
left_spacing=80; % pixel, the space between subplot and the left things
right_spacing=50; % pixel, the space right of the last column of subplot
upper_spacing=70; % pixel, the space between subplot and the upper things
lower_spacing=30; % pixel, the space below the legend box
legend_height=100; % pixel, the height of legend box

plot_n_col=3; % the column number of subplot
plot_n_row=2; % the row number of subplot

%% main
subject_error_arr=[];

mkdir(input_dir,'arrangement');
for sbj=1:length(target_name_arr)
    % make output folder and find fitting folder
    mkdir(fullfile(input_dir,'arrangement',target_name_arr{sbj}));
    temp_target_folder=dir(fullfile(input_dir,[target_name_arr{sbj} '*']));
    assert(length(temp_target_folder)==1,['Error! more than 1 fitting folder for ' target_name_arr{sbj} '!']);
    temp_target_folder=fullfile(input_dir,temp_target_folder.name);
    
    % load ANN model
    ANN_model=load(fullfile(model_dir,[model_name_arr{sbj} '_model.mat'])); % net, param_range
    net=ANN_model.net; param_range=ANN_model.param_range;
    
    target_spec=load(fullfile(target_dir,[target_name_arr{sbj} '.txt']));
    fitting_wl=load(fullfile(temp_target_folder,'fitting_1','wl.txt'));
    init_param_arr=load(fullfile(temp_target_folder,'init_arr.txt'));
    fitting_result=[];
    fitting_op_arr=[];
    
    lambda=fitting_wl;
    fun_init_param_to_mu_spec(); % load the epsilon for fitting wl
    
    interped_target_spec=interp1(target_spec(:,1),target_spec(:,2:end),fitting_wl); % for calculate error
    
    fitted_ann_spec_arr=[];
    init_ann_spec_arr=[];
    
    for j=1:times_to_fitting
%         if j==1
%             SDS_choosed=load(fullfile(temp_target_folder,['fitting_' num2str(j)],'setting_record.mat'),'SDS_choosed');
%             SDS_choosed=SDS_choosed.SDS_choosed; % the SDS chosen to fitted in the target spectrum
%         end
        
        fprintf('Process %s fitting %d\n',target_name_arr{sbj},j);
        
        % load the fitting result
        temp_route=load(fullfile(temp_target_folder,['fitting_' num2str(j)],'fitRoute.txt'));
        fitting_result(j,:)=temp_route(end,:);
        fitting_op_arr(:,:,j)=load(fullfile(temp_target_folder,['fitting_' num2str(j)],'fitted_mu.txt'));

        % forward the fitted spec and the init spec
        [OP_arr,~]=fun_param_to_mu(fitting_result(j,1:num_fitted_param),0);
        fitted_ann_spec=fun_ANN_forward(OP_arr);
        fitted_ann_spec_arr(:,:,j)=fitted_ann_spec(:,SDS_sim_correspond);
        
        % save the fitted OPs
        to_save=[lambda OP_arr];
        save(fullfile(input_dir,'arrangement',target_name_arr{sbj},['fitted_OP_' num2str(j) '.txt']),'to_save','-ascii','-tabs');
        to_save=[lambda fitted_ann_spec];
        save(fullfile(input_dir,'arrangement',target_name_arr{sbj},['fitted_spec_' num2str(j) '.txt']),'to_save','-ascii','-tabs');
        
        [OP_arr,~]=fun_param_to_mu(init_param_arr(j,1:num_fitted_param),0);
        init_ann_spec=fun_ANN_forward(OP_arr);
        init_ann_spec_arr(:,:,j)=init_ann_spec(:,SDS_sim_correspond);
        
        % calculate the error
        calculated_init_error=sqrt(mean((init_ann_spec_arr(:,:,j)./interped_target_spec-1).^2,1));
        calculated_fitted_error=sqrt(mean((fitted_ann_spec_arr(:,:,j)./interped_target_spec-1).^2,1));
        fitting_result(j,num_fitted_param+1:num_fitted_param+num_SDS+1)=[calculated_fitted_error sqrt(mean(calculated_fitted_error(SDS_choosed).^2))];
        init_param_arr(j,num_fitted_param+1:num_fitted_param+num_SDS+1)=[calculated_init_error sqrt(mean(calculated_init_error(SDS_choosed).^2))];
    end

    %% output CSV
    [~,error_index]=sort(fitting_result(:,end),'ascend');
    fitting_result_sort=[error_index fitting_result(error_index,:)];
    init_param_arr_sort=init_param_arr(error_index,:);
    save(fullfile(input_dir,'arrangement',[target_name_arr{sbj} '_fitting_result_sort.txt']),'fitting_result_sort','-ascii','-tabs');
    save(fullfile(input_dir,'arrangement',[target_name_arr{sbj} '_init_param_arr_sort.txt']),'init_param_arr_sort','-ascii','-tabs');

    subject_error_arr(1,sbj)=fitting_result_sort(1,end);
    
    fid=fopen(fullfile(input_dir,'arrangement',[target_name_arr{sbj} '_fitting_results.csv']),'w');
    % header
    fprintf(fid,',A1,K1,A2,K2,A4,K4,hc1,sto2_1,hc2,sto2_2,hc4,sto2_4,mel');
    for s=1:num_SDS
        fprintf(fid,',error %d',s);
    end
    fprintf(fid,',error_choosed,,A1,K1,A2,K2,A4,K4,hc1,sto2_1,hc2,sto2_2,hc4,sto2_4,mel');
    for s=1:num_SDS
        fprintf(fid,',error %d',s);
    end
    fprintf(fid,',error_choosed\n');
    for j=error_index'
        fprintf(fid,'Fitting %d',j);
        fprintf(fid,',%f',fitting_result(j,1:num_fitted_param));
        fprintf(fid,',%.2f%%',fitting_result(j,num_fitted_param+1:end)*100);
        fprintf(fid,',init param');
        fprintf(fid,',%f',init_param_arr(j,1:num_fitted_param));
        fprintf(fid,',%.2f%%',init_param_arr(j,num_fitted_param+1:end)*100);
        fprintf(fid,'\n');
    end
    fclose(fid);
    
    if do_plot_anyPlot
        %% plot the OPs
        fprintf('Plot fitted OP\n');
        legend_arr={};
        for j=1:times_to_fitting
            legend_arr{j}=['rank ' num2str(times_to_fitting+1-j)];
        end

        fig=figure('Units','pixels','position',[0 0 (left_spacing+subplot_width)*plot_n_col+right_spacing (upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing+lower_spacing]);
        set(fig,'visible','off');
        op_name={'\mu_{a','\mu_{s'};
        layer_name_arr={'scalp','skull','CSF','GM'};
        colormap_arr=jet(times_to_fitting);
        colormap_arr=colormap_arr(end:-1:1,:);
        subplot_index=1;
        for L=[1 2 4]
            for opi=1:2
%                 row_index=ceil(subplot_index/plot_n_col);
%                 col_index=subplot_index-(row_index-1)*plot_n_col;
                row_index=opi;
                col_index=L;
                if L==4
                    col_index=3;
                end
                axes;
                subplot_index=subplot_index+1;
                hold on;
                for j=times_to_fitting:-1:1
                    plot(fitting_wl,fitting_op_arr(:,2*(L-1)+opi,error_index(j)),'Color',colormap_arr(j,:),'LineWidth',lineWidth);
                end
                % scale to ub and lb
                ylim([param_range(2,(opi-1)*4+L) param_range(1,(opi-1)*4+L)]);
                title([op_name{opi} ',' layer_name_arr{L} '}']);
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
        text(((left_spacing+subplot_width)*plot_n_col+right_spacing)/2,(upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing*2/3+lower_spacing,strrep(target_name_arr{sbj},'_',' '),'Unit','pixels','fontsize',fontSize, 'FontName', 'Times New Roman','HorizontalAlignment','center','VerticalAlignment','middle')

        print(fullfile(input_dir,'arrangement',[target_name_arr{sbj} '_mu.png']),'-dpng','-r200');
        close all;
    
        %% plot the fitted spec together
        fprintf('Plot the fitted spec\n');
        legend_arr={'target'};
        for j=1:times_to_fitting
            legend_arr{end+1}=['rank ' num2str(times_to_fitting+1-j) ', ' num2str(fitting_result(error_index(times_to_fitting+1-j),end)*100,'%.2f%%')];
        end

        fig=figure('Units','pixels','position',[0 0 (left_spacing+subplot_width)*plot_n_col+right_spacing (upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing+lower_spacing]);
        set(fig,'visible','off');
        colormap_arr=jet(times_to_fitting);
        colormap_arr=colormap_arr(end:-1:1,:);
        subplot_index=1;
        for s=1:num_SDS
            row_index=ceil(subplot_index/plot_n_col);
            col_index=subplot_index-(row_index-1)*plot_n_col;
            subplot(2,3,subplot_index);
            subplot_index=subplot_index+1;
            hold on;
            plot(target_spec(:,1),target_spec(:,s+1),'LineWidth',lineWidth);
            for j=times_to_fitting:-1:1
                plot(lambda,fitted_ann_spec_arr(:,s,error_index(j)),'--','Color',colormap_arr(j,:),'LineWidth',lineWidth);
            end
            title(['SDS = ' num2str(SDS_dist_arr(s)) ' cm']);
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
        axes;
        axis off;
        set(gca,'Unit','pixels','Position',[0 0 (left_spacing+subplot_width)*plot_n_col+right_spacing (upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing+lower_spacing]);
        text(((left_spacing+subplot_width)*plot_n_col+right_spacing)/2,(upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing*2/3+lower_spacing,strrep(target_name_arr{sbj},'_',' '),'Unit','pixels','fontsize',fontSize, 'FontName', 'Times New Roman','HorizontalAlignment','center','VerticalAlignment','middle')

        print(fullfile(input_dir,'arrangement',[target_name_arr{sbj} '_spec.png']),'-dpng','-r200');
        close all;
    end
    
    %% plot the fitted spec
    if do_plot_individual_fitting && do_plot_anyPlot
%         for jj=1:times_to_fitting % the rank
        for jj=1:1 % the rank
            fprintf('Plot %s fitting %d\n',target_name_arr{sbj},jj);
            j=error_index(jj); % the index of fitting
            fig=figure('Units','pixels','position',[0 0 1600 900]);
            set(fig,'visible','off');
            ti=tiledlayout('flow','TileSpacing','compact','Padding','none');
            for s=1:num_SDS
                nexttile();
                hold on;
                plot(target_spec(:,1),target_spec(:,s+1),'LineWidth',lineWidth);
                plot(lambda,init_ann_spec_arr(:,s,j),'LineWidth',lineWidth);
                plot(lambda,fitted_ann_spec_arr(:,s,j),'LineWidth',lineWidth);
                grid on;

                lgd=legend({'target',['init, err=' num2str(init_param_arr(j,num_fitted_param+s)*100,'%.2f%%')],['fitted, err=' num2str(fitting_result(j,num_fitted_param+s)*100,'%.2f%%')]},'Location','southoutside');
                lgd.NumColumns=3;
                title(['SDS ' num2str(s)]);
                set(gca,'fontsize',fontSize, 'FontName', 'Times New Roman');
            end
            title(ti,[strrep(target_name_arr{sbj},'_',' ') ' rank ' num2str(jj) ' = fitting ' num2str(j) ', fitting error = ' num2str(fitting_result(j,end)*100,'%.2f%%')], 'FontName', 'Times New Roman');
            print(fullfile(input_dir,'arrangement',target_name_arr{sbj},['fitted_spec_r' num2str(jj) '_f' num2str(j) '.png']),'-dpng','-r200');
            close all;
        end
    end
end

save(fullfile(input_dir,'arrangement','subject_error_arr.txt'),'subject_error_arr','-ascii','-tabs');

disp('Done!');