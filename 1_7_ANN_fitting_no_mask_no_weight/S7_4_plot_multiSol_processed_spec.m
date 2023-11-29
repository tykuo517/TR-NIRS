%{
Plot the subject fitted spec after process the multiple solution

Benjamin Kao
Last update: 2021/03/31
%}

clc;clear;close all; clearvars -global;

%% param
input_dir='fitted_result';
target_name_arr={'tc_mean','ww_mean','wh_mean','yf_mean','kb_mean'}; % the name of the target spectrum
to_process_sbj_index=3; % the index of the subject to process
target_dir='input_target_2'; % the folder containing the target spectrum

fitting_dir={'fitting_SDS1234_3','fitting_SDS2345_3','fitting_SDS123456_3','fitting_SDS12345_3','fitting_SDS346_3'}; % please move the fitting folders into this folder first.
sbj_fitting_dir_arr=[3 5 1 4 1]; % which fitting folder does this subject use
num_fitted_param=13; % the number of the fitted parameters

num_SDS=6;

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

%% main
sbj_name=target_name_arr{to_process_sbj_index};
choosed_rank=load(fullfile(input_dir,[sbj_name '_choosed_rank.txt']));

fitted_spec_arr=[];
for i=1:length(choosed_rank)
    fitted_spec_arr(:,:,i)=load(fullfile(input_dir,[sbj_name '_fitted_spec_' num2str(i) '.txt']));
end

fitting_result=load(fullfile(fitting_dir{sbj_fitting_dir_arr(to_process_sbj_index)},'arrangement',[sbj_name '_fitting_result_sort.txt']));
fitting_result=fitting_result(choosed_rank,2:end);

%% plot the fitted spec together
target_spec=load(fullfile(target_dir,[target_name_arr{to_process_sbj_index} '.txt']));

fprintf('Plot the fitted spec\n');

fig=figure('Units','pixels','position',[0 0 (left_spacing+subplot_width)*plot_n_col+right_spacing (upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing+lower_spacing]);
set(fig,'visible','off');
colormap_arr=lines(length(choosed_rank)+1);

SDS_cm_arr=[0.8 1.5 2.12 3 3.35 4.5];
ti=tiledlayout('flow','TileSpacing','compact','Padding','none');
for s=1:num_SDS
    nexttile();
    hold on;
    plot(target_spec(:,1),target_spec(:,s+1),'Color',colormap_arr(1,:),'LineWidth',lineWidth);
    for j=1:length(choosed_rank)
        plot(fitted_spec_arr(:,1,j),fitted_spec_arr(:,s+1,j),'--','Color',colormap_arr(j+1,:),'LineWidth',lineWidth);
    end
    xlabel('wavelength(nm)');
    ylabel('reflectance');
    title(['SDS = ' num2str(SDS_cm_arr(s)) ' cm']);
    set(gca,'fontsize',fontSize, 'FontName', 'Times New Roman');
    grid on;
    yylim=ylim(); yylim(1)=0; ylim(yylim);

    legend_arr={'target'};
    for j=1:length(choosed_rank)
        legend_arr{end+1}=['rank ' num2str(choosed_rank(j)) ', spec error= ' num2str(fitting_result(j,num_fitted_param+s)*100,'%.2f%%')]; % spec error
    end
    lgd=legend(legend_arr,'Location','southoutside','fontsize',lgdFontSize);
    lgd.NumColumns=lgdNumCol;
end

print(fullfile(input_dir,[target_name_arr{to_process_sbj_index} '_choosed_spec.png']),'-dpng','-r200');
close all;

disp('Done!');