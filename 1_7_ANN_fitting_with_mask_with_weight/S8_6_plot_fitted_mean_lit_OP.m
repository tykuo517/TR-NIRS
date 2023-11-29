%{
Plot the fitted OP, the mean of fitted Op and the literature OP

Benjamin Kao
Last update: 2021/03/25
%}

clc;clear;close all;

global lambda net param_range 

%% param
input_dir='fitted_result';
output_dir='saved_OPs';

litOP_dir=fullfile('literature_OPs','OPs_to_sim_15');

target_name_arr={'tc_mean','ww_mean','wh_mean','yf_mean','kb_mean'}; % the name of the target spectrum

subplot_height=300; % pixel, the height of subplot
subplot_width=450; % pixel, the width of subplot
left_spacing=100; % pixel, the space between subplot and the left things
right_spacing=50; % pixel, the space right of the last column of subplot
upper_spacing=110; % pixel, the space between subplot and the upper things
lower_spacing=30; % pixel, the space below the legend box
legend_height=50; % pixel, the height of legend box

plot_n_col=3; % the column number of subplot
plot_n_row=2; % the row number of subplot

lineWidth=2;
fontSize=18;
lgdFontSize=18;
lgdNumCol=11;

%% main

%% load the subject OP
all_sbj_OP_arr=[];
all_answer_index_arr=[];
all_answer_subject_arr=[];

for i=1:length(target_name_arr)
    sbj_answer_number=load(fullfile(input_dir,output_dir,[target_name_arr{i} '_OP_info.mat']),'sbj_answer_number');
    sbj_answer_number=sbj_answer_number.sbj_answer_number;
    
    % load the fitted OP
    for j=1:sbj_answer_number
        fitted_OP_arr=load(fullfile(input_dir,output_dir,[target_name_arr{i} '_OP_' num2str(j) '.txt']));
        if length(all_sbj_OP_arr)==0
            all_sbj_OP_arr(:,:,1)=fitted_OP_arr;
        else
            all_sbj_OP_arr(:,:,end+1)=fitted_OP_arr;
        end
        all_answer_subject_arr(end+1,1)=i;
        all_answer_index_arr(end+1,1)=j;
    end
end

%% load the mean fitted OP
mean_OP=load(fullfile(input_dir,output_dir,[target_name_arr{1} '_OP_' num2str(sbj_answer_number+1) '.txt']));

%% load the literature OP
literature_wl=load(fullfile(litOP_dir,'sim_wl.txt'));
literature_op=load(fullfile(litOP_dir,'toSim_OP_1.txt'));
litOP=[literature_wl literature_op];

%% plot

legend_arr={'average fitted OP','average literature OP'};
for i=1:size(all_sbj_OP_arr,3)
    if all_answer_index_arr(i)==1
        legend_arr{end+1}=['subject ' num2str(all_answer_subject_arr(i))];
    else
        legend_arr{end+1}=['subject ' num2str(all_answer_subject_arr(i)) '_' num2str(all_answer_index_arr(i))];
    end
end

fig=figure('Units','pixels','position',[0 0 (left_spacing+subplot_width)*plot_n_col+right_spacing (upper_spacing+subplot_height)*plot_n_row+legend_height+upper_spacing+lower_spacing]);
set(fig, 'visible', 'off');
OP_name_arr={'\mu_{a,scalp}','\mu_{s,scalp}','\mu_{a,skull}','\mu_{s,skull}','\mu_{a,CSF}','\mu_{s,CSF}','\mu_{a,GM}','\mu_{s,GM}',};

s=1;
for L=[1 2 3 4 7 8]
    row_index=ceil(s/plot_n_col);
    col_index=s-(row_index-1)*plot_n_col;
    axes;
    
    hold on;
    plot(mean_OP(:,1),mean_OP(:,L+1),'LineWidth',3);
    plot(litOP(:,1),litOP(:,L+1),'LineWidth',3);
    for i=1:size(all_sbj_OP_arr,3)
        plot(all_sbj_OP_arr(:,1,i),all_sbj_OP_arr(:,L+1,i),'LineWidth',1.5);
    end
    title(OP_name_arr{L});
    xlabel('wavelength(nm)');
    ylabel([OP_name_arr{L} '(1/cm)']);
    grid on;
    xlim([700 900])
    
    if s==6
        lgd=legend(legend_arr,'Location','southoutside','fontsize',lgdFontSize);
        lgd.NumColumns=lgdNumCol;
        set(lgd,'Unit','pixels','position',[left_spacing lower_spacing plot_n_col*subplot_width+(plot_n_col-1)*left_spacing legend_height]);
    end
    
    set(gca,'fontsize',fontSize, 'FontName', 'Times New Roman');
    set(gca,'Unit','pixels','Position',[left_spacing+(left_spacing+subplot_width)*(col_index-1) lower_spacing+legend_height+upper_spacing+(subplot_height+upper_spacing)*(plot_n_row-row_index) subplot_width subplot_height]);
    s=s+1;
end

print(fullfile(input_dir,output_dir,'mean_OP.png'),'-dpng','-r200');

disp('Done!');