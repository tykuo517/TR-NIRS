%{
Calculate mean and std OP error of each fitting result and plot in bar
chart

Ting-Yi Kuo
Last update: 2023/05/11
%}

clc;clear;%close all;

global lambda Lbound Ubound net param_range;

%% param

Add_error_mode=1; % if =0, use result without error; =1 to use results with error

input_dir='test_fitting_2023-07-17-16-40-04_MC_5'; % please move the fitting folders into this folder first.
subject_name_arr={'KB'};
num_anser_to_generate=10; % number of target spec (true answer)
num_error_to_generate=15; % number of adding noise to the same, the first one will have no error
times_to_fitting=20; % number of fitting using different init value

fitting_dir={'fitting_SDS345_gate1-5'}; % the fitting folder
num_SDS=5; % how many SDS are in the target spectrum
num_gate=10;
num_fitted_param=5; % the number of the fitted parameters

OP_total_result=zeros(2,4,length(fitting_dir));

if Add_error_mode==0
    to_process_fitting_index=1;
elseif Add_error_mode==1
    to_process_fitting_index=2:num_error_to_generate;
end

OP_result=[];

for sbj=1:length(subject_name_arr)
    for dir=1:length(fitting_dir)
        % make output folder and find fitting folder
        mkdir(fullfile(input_dir,'arrangement',subject_name_arr{sbj},fitting_dir{dir}));
        fitting_op_arrange=[];

        for target_i=1:num_anser_to_generate
            for error_i=to_process_fitting_index

                temp_target_folder=fullfile(input_dir,subject_name_arr{sbj},['target_' num2str(target_i) '_' num2str(error_i)],fitting_dir{dir});

                fitting_op_arrange_temp=load(fullfile(temp_target_folder,'fitting_result_sort.txt'));
                fitting_op_arrange(end+1,:)=fitting_op_arrange_temp(1,num_fitted_param+num_SDS*num_gate+3:end); %num_fitted_param+num_SDS+3:end
            end
        end
        OP_result(1:sbj*num_anser_to_generate*length(to_process_fitting_index),1:4,dir)=[OP_result(:,:,dir); fitting_op_arrange];
%         OP_total_result(:,:,dir)=OP_result;
    %     save(fullfile(input_dir,'arrangement',subject_name_arr{sbj},fitting_dir,'OP_reult_mean_std.txt'),'OP_result','-ascii','-tabs');
    %     save(fullfile(input_dir,'arrangement',subject_name_arr{sbj},fitting_dir,'OP_reult.txt'),'fitting_op_arrange','-ascii','-tabs');
    end
end

OP_mean_result=zeros(length(fitting_dir),4);
OP_mean_result=mean(OP_result,1);
OP_std_result=std(OP_result(:,:,dir),[],1);

%% Plot
fitting_choosed={};
for dir=1:length(fitting_dir)
    dir_part=strsplit(fitting_dir{dir},'_');
    fitting_choosed{end+1}=[dir_part{2} ' ' dir_part{3}];
end

colormap_arr=jet(length(fitting_dir));
colormap_arr=colormap_arr(end:-1:1,:);
    
figure('Position',[0 0 840 680]);
y=100*OP_mean_result';
b=bar(y,'grouped');
hold on;
[ngroups,nbars] = size(y);
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',y,100*OP_std_result','k','linestyle','none');

set(gca, 'XTickLabel', {'\mu_{a,skull}' '\mu_{a,GM}' '\mu_{s,GM}' 'all OPs'});
ylabel('error (%)');
title('OP errors of different combinations of SDS and gate')
legend(fitting_choosed,'Location','northwest');

if Add_error_mode==0
    print(fullfile(input_dir,'arrangement','opErrorResult_noError_bar.png'),'-dpng','-r200');
elseif Add_error_mode==1
    print(fullfile(input_dir,'arrangement','opErrorResult_withError_bar.png'),'-dpng','-r200');
end



disp('Done!');
            
            