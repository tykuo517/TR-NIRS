%{
Plot the scattering plot of the fitted tissue parameter and the spectral error

Benjamin Kao
Last update: 2021/01/25
%}

clc;clear;close all;

%% param
input_dir='fitting_SDS12345_3'; % please move the fitting folders into this folder first.

target_name_arr={'tc_mean','ww_mean','wh_mean','yf_mean','kb_mean'}; % the name of the target spectrum
model_name_arr={'ZJ','WW2','WH2','YF','KB'}; % the name of ANN model corresponding to each target spec

num_SDS=6; % how many SDS are in the target spectrum
SDS_sim_correspond=[1 2 3 4 5 6]; % the SDS index in the simulated spectrum corresponding to each SDS in the target spectrum, SDS_sim_correspond(x)=y means 'measure x = sim y'

param_name_arr={'hc1','sto2_1','hc2','sto2_2','hc4','sto2_4','mel'};
num_fitted_param=13; % the number of the fitted parameters
times_to_fitting=20; % number of init value used to fitting

model_dir='model_arrange'; % the folder containing the arranged ANN file
target_dir='input_target_2'; % the folder containing the target spectrum

fontSize=16;
lineWidth=2;
lgdFontSize=14;
lgdNumCol=10;

subplot_height=330; % pixel, the height of subplot
subplot_width=450; % pixel, the width of subplot
left_spacing=100; % pixel, the space between subplot and the left things
right_spacing=50; % pixel, the space right of the last column of subplot
upper_spacing=100; % pixel, the space between subplot and the upper things
lower_spacing=0; % pixel, the space below the legend box
legend_height=0; % pixel, the height of legend box

plot_n_col=4; % the column number of subplot
plot_n_row=2; % the row number of subplot

%% init
mua_coef_bound=load(fullfile(input_dir,'mua_coef_bound.txt'));

%% main
for sbj=1:length(target_name_arr)
    fprintf('Plot %s\n',target_name_arr{sbj});
    
    fitted_param_arr=load(fullfile(input_dir,'arrangement',[target_name_arr{sbj} '_fitting_result_sort.txt']));
    fitted_param_arr=fitted_param_arr(:,2:end); % delete the 1st column of fitting index
    init_param_arr=load(fullfile(input_dir,'arrangement',[target_name_arr{sbj} '_init_param_arr_sort.txt']));
    
    % plot the relationship between each parameter and error
    fig=figure('Units','pixels','position',[0 0 (left_spacing+subplot_width)*plot_n_col+right_spacing (upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing+lower_spacing]);
    set(fig,'visible','off');
    subplot_index=1;
    for p=1:7
        row_index=ceil(subplot_index/plot_n_col);
        col_index=subplot_index-(row_index-1)*plot_n_col;
        subplot(plot_n_row,plot_n_col,subplot_index);
        subplot_index=subplot_index+1;
        plot(fitted_param_arr(:,end),fitted_param_arr(:,6+p),'o','LineWidth',lineWidth);
        % scale to ub and lb
        ylim([mua_coef_bound(1,p) mua_coef_bound(2,p)]);
        
        xlabel('error');
        
        title(param_name_arr{p});
        set(gca,'fontsize',fontSize, 'FontName', 'Times New Roman');
        set(gca,'Unit','pixels','Position',[left_spacing+(left_spacing+subplot_width)*(col_index-1) lower_spacing+legend_height+upper_spacing+(subplot_height+upper_spacing)*(plot_n_row-row_index) subplot_width subplot_height]);
        grid on;
    end
    
    % plot the relationship between init error and fitted error
    row_index=ceil(subplot_index/plot_n_col);
    col_index=subplot_index-(row_index-1)*plot_n_col;
    subplot(plot_n_row,plot_n_col,subplot_index);
    subplot_index=subplot_index+1;
    plot(init_param_arr(:,end),fitted_param_arr(:,end),'o','LineWidth',lineWidth);
    xlabel('init error');
    ylabel('fitted error');
    set(gca,'fontsize',fontSize, 'FontName', 'Times New Roman');
    set(gca,'Unit','pixels','Position',[left_spacing+(left_spacing+subplot_width)*(col_index-1) lower_spacing+legend_height+upper_spacing+(subplot_height+upper_spacing)*(plot_n_row-row_index) subplot_width subplot_height]);
    grid on;
    
    % add the title
    axes;
    axis off;
    set(gca,'Unit','pixels','Position',[0 0 (left_spacing+subplot_width)*plot_n_col+right_spacing (upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing+lower_spacing]);
    text(((left_spacing+subplot_width)*plot_n_col+right_spacing)/2,(upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing*2/3+lower_spacing,strrep(target_name_arr{sbj},'_',' '),'Unit','pixels','fontsize',fontSize, 'FontName', 'Times New Roman','HorizontalAlignment','center','VerticalAlignment','middle')
    
    print(fullfile(input_dir,'arrangement',[target_name_arr{sbj} '_paramError.png']),'-dpng','-r200');
    close all;
end

disp('Done!');