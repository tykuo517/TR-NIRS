%{
arrange the calculated OP error meam, std or confidence interval

Benjamin Kao
Last update: 2021/03/29
%}

clc;clear;close all; clearvars -global;

%% param
fitting_dir_arr={'test_fitting_2023-10-31-21-57-22'};
fitting_SDS_dir_arr={'fitting_SDS123456_345_gate1-5','fitting_SDS123456__gate1-5'};
fitting_dir_index=[1 1];
add_noise_arr=[1 1];
fitting_name_arr={'CW+TR','CW'}; %'SDS 12 merged no noise',

%% main
std_arr=[];
mean_arr=[];
RMSE_arr=[];
ci_arr={};
for i=1:3
    ci_arr{i}=[];
end

for i=1:length(fitting_SDS_dir_arr)
    if add_noise_arr(i)==0
        to_load_name='OP_error_arr_noError.mat';
    else
        to_load_name='OP_error_arr_Error.mat';
    end
    error_info=load(fullfile(fitting_dir_arr{fitting_dir_index(i)},'arrangement',fitting_SDS_dir_arr{i},to_load_name));
    
    RMSE_arr(i,:)=sqrt(mean(error_info.OP_error_arr.^2,[1 3]));
    
    for opi=1:6
        for ci_l=1:3 % confidence level
            for j=1:2
                ci_arr{ci_l}(i,opi,j)=error_info.OP_error_CI{opi}(ci_l,j);
            end
        end
        std_arr(i,opi)=error_info.total_std_error_arr(opi);
        mean_arr(i,opi)=error_info.total_mean_OP_error_arr(opi);
    end
end

%% output the mean error
fprintf('mean error:\n');
fprintf('| op | mua1 | mus1 | mua2 | mus2 | mua4 | mus4 |\n|:------ | ---- | ---- | ---- |:---- | ---- |:---- |\n'); % header
for i=1:length(fitting_SDS_dir_arr)
    fprintf('|%s',fitting_name_arr{i});
    for opi=1:6
        fprintf('|%.2f%%',mean_arr(i,opi)*100);
    end
    fprintf('|\n');
end

%% output the error std
fprintf('\nerror std:\n');
fprintf('| op | mua1 | mus1 | mua2 | mus2 | mua4 | mus4 |\n|:------ | ---- | ---- | ---- |:---- | ---- |:---- |\n'); % header
for i=1:length(fitting_SDS_dir_arr)
    fprintf('|%s',fitting_name_arr{i});
    for opi=1:6
        fprintf('|%.2f%%',std_arr(i,opi)*100);
    end
    fprintf('|\n');
end

%% output the RMSPE
fprintf('\nRMSPE:\n');
fprintf('| op | mua1 | mus1 | mua2 | mus2 | mua4 | mus4 |\n|:------ | ---- | ---- | ---- |:---- | ---- |:---- |\n'); % header
for i=1:length(fitting_SDS_dir_arr)
    fprintf('|%s',fitting_name_arr{i});
    for opi=1:6
        fprintf('|%.2f%%',RMSE_arr(i,opi)*100);
    end
    fprintf('|\n');
end

% plot RMSE in bar chart
% RMSE_arr=[RMSE_arr; 0.1839 0.0747 0.3852 0.1419 0.2665 0.6902];
figure('Position',[0 0 840 680]);
y=100*RMSE_arr';
b=bar(y,'grouped');
xticklabels({'\mu_{a,scalp}','\mu_{s,scalp}','\mu_{a,skull}','\mu_{s,skull}','\mu_{a,GM}','\mu_{s,GM}'});
ylabel('RMSPE (%)');
legend([fitting_name_arr 'CW'],'Location','northwest');
mkdir('results')
print(fullfile('results','OP_result_RMSE_bar.png'),'-dpng','-r200');

%% output the confidence interval
fprintf('\nconfidence interval:\n');
for ci_l=1:3
    fprintf('confidence level %f:\n',error_info.confidence_arr(ci_l));
    fprintf('| op | mua1 | mus1 | mua2 | mus2 | mua4 | mus4 |\n|:------ | ---- | ---- | ---- |:---- | ---- |:---- |\n'); % header
    for i=1:length(fitting_SDS_dir_arr)
        fprintf('|%s',fitting_name_arr{i});
        for opi=1:6
            fprintf('|');
            for j=1:2
                fprintf('%.2f%%',ci_arr{ci_l}(i,opi,j)*100);
                if j==1
                    fprintf('~');
                end
            end
        end
        fprintf('|\n');
    end
end