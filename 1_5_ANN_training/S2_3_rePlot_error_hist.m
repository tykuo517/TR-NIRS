%{
Replot the training error histogram

Benjamin Kao
Last update: 2020/12/23
%}

clc;clear;close all;

%% param
subject_name_arr={'KB'}; % the name of the subject
input_dir_arr={'KB_2023-10-11-20-35-47'};
num_SDS=5;
num_gate=10;
SDS_dist_arr=[1.5 2.2 2.9 3.6 4.3]; % cm
fontSize=10;
marker_size=20;
line_width=2;

for sbj_i=1:length(input_dir_arr)
    input_dir=input_dir_arr{sbj_i};
    %% init
    fprintf('Processing %s\n',input_dir);
    load(fullfile(input_dir,'ANN__train_info.mat'));
    
    error_mean=mean(abs(error),1);
    error_std=std(error,[],1);
    error_rmspe=sqrt(mean(error.^2,1));
    
    %% main
    for s=1:num_SDS
        fig=figure('Units','pixels','position',[0 0 1980 1080]);
        set(fig, 'visible', 'off');
        ti=tiledlayout('flow','TileSpacing','compact','Padding','none');
        for g=1:num_gate
            nexttile();
            [a,b]=hist(error(:,num_gate*(s-1)+g),50);
            bar(b,a+1);
            set(gca, 'YScale', 'log');
            grid on;
            xlabel('Error');
            ylabel('testing data number');
            title({['Gate = ' num2str(g)],['mean=' num2str(error_mean(num_gate*(s-1)+g)*100,'%.2f%%') ', std=' num2str(error_std(num_gate*(s-1)+g)*100,'%.2f%%') ', rmspe=' num2str(error_rmspe(num_gate*(s-1)+g)*100,'%.2f%%')]});
            set(gca,'fontsize',fontSize, 'FontName', 'Times New Roman');
            fprintf('---- error= %2.4f%%, rmse= %2.4f%%, std=%2.4f%%\n',error_mean(num_gate*(s-1)+g)*100,error_rmspe(num_gate*(s-1)+g)*100,error_std(num_gate*(s-1)+g)*100);
        end
        title(ti,[subject_name_arr{sbj_i} '  SDS ' num2str(SDS_dist_arr(s)) ' cm training result'], 'FontName', 'Times New Roman');

        print(fullfile(input_dir,['SDS_' num2str(SDS_dist_arr(s)) 'cm_testing_error_histogram.png']),'-dpng','-r200');
    end
    close all;
end

disp('Done!');