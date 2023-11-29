%{
Calculate the fitted OP error

Benjamin Kao
Last update: 2021/03/16
%}

clc;clear;close all; clearvars -global;

global lambda net param_range 

%% param
error_range_threshold=0.005; % plot the fitting in this distance to the best result
max_choose_number=1; % the max number of fitting

input_dir='test_fitting_2023-10-31-21-57-22'; % please move the fitting folders into this folder first.
subject_name_arr={'KB'};%'KB','WH','WW'
do_use_add_error=1; % if =1, calculate the OP error of the added noise results; if =0, calculate the error without noise
num_anser_to_generate=10; % number of target spec (true answer)
num_error_to_generate=15; % number of adding noise to the same, the first one will have no error
times_to_fitting=20; % number of fitting using different init value

fitting_dir='fitting_SDS123456_345_gate1-5';

% CW setting
num_SDS_cw=6; % how many SDS are in the target spectrum
SDS_choosed_cw=[1 2 3 4 5 6]; % the SDS chosen to fitted in the target spectrum
SDS_sim_correspond_cw=[1 2 3 4 5 6]; % the SDS index in the simulated spectrum corresponding to each SDS in the target spectrum, SDS_sim_correspond(x)=y means 'measure x = sim y'

% TR setting
num_SDS_tr=5;
SDS_choosed_tr=[]; % the SDS chosen to fitted in the target spectrum
SDS_sim_correspond_tr=[1 2 3 4 5]; % the SDS index in the simulated spectrum corresponding to each SDS in the target spectrum, SDS_sim_correspond(x)=y means 'measure x = sim y'

num_gate=10;
gate_choosed=[1 2 3 4 5];
gate_sim_correspond=[1 2 3 4 5 6 7 8 9 10];

num_fitted_param=13; % the number of the fitted parameters

fitting_wl_file='fitting_wl.txt'; % the file in the epsilon folder containing the fitting wavelength
fitting_wl_tr=810;

model_dir='model_arrange'; % the folder containing the arranged ANN file

do_plot_anyPlot=0; % =1 to plot any plot
do_plot_individual_fitting=1; % =1 to plot the result of each individual fitting

fontSize=14;
lineWidth=2;
lgdFontSize=14;
lgdNumCol=10;

subplot_height=300; % pixel, the height of subplot
subplot_width=450; % pixel, the width of subplot
left_spacing=120; % pixel, the space between subplot and the left things
right_spacing=50; % pixel, the space right of the last column of subplot
upper_spacing=150; % pixel, the space between subplot and the upper things
lower_spacing=30; % pixel, the space below the legend box
legend_height=0; % pixel, the height of legend box

plot_n_col=7; % the column number of subplot
plot_n_row=4; % the row number of subplot


%% init
fitting_wl=load(fullfile('epsilon',fitting_wl_file));
fitting_wl=[fitting_wl; fitting_wl_tr];
unique(fitting_wl);
lambda=fitting_wl;

mkdir(fullfile(input_dir,'arrangement',fitting_dir));

OP_column_arr=[1 2 3 4 7 8]; % the column in OP array of each OP

% load the param and OP answer
% param_answer_arr=load(fullfile(input_dir,'answers','param_answer_arr.txt'));
OP_answer_arr=[];
for i=1:num_anser_to_generate
    temp_OP_arr=load(fullfile(input_dir,'answers',['OP_ans_' num2str(i) '.txt']));
    OP_answer_arr(:,:,i)=interp1(temp_OP_arr(:,1),temp_OP_arr(:,1+OP_column_arr),lambda);
end

colormap_arr=jet(num_error_to_generate);

legend_arr={};
for error_i=1:num_error_to_generate
    legend_arr{error_i}=['error ' num2str(error_i)];
end

%% main
OP_error_arr=[]; % wl * OP * target
OP_error_counter=1;

for sbj=1:length(subject_name_arr)
% for sbj=[4]
    % make output folder and find fitting folder
    

    % load the fitting result
    for target_i=1:num_anser_to_generate
        fprintf('Process %s target %d\n',subject_name_arr{sbj},target_i);
        
        fitted_result_arr=[];
        
        if do_use_add_error
            to_process_fitting_index=2:num_error_to_generate;
        else
            to_process_fitting_index=1;
        end
        
        for error_i=to_process_fitting_index

            temp_target_folder=fullfile(input_dir,subject_name_arr{sbj},['target_' num2str(target_i) '_' num2str(error_i)],fitting_dir);
            
            temp_fitting_result=load(fullfile(temp_target_folder,'fitting_result_sort.txt'));
            orig_init_index=temp_fitting_result(:,1);
            temp_fitting_result=temp_fitting_result(:,2:end);
            fitted_result_arr(:,:,error_i)=temp_fitting_result;
            
            % choose the closest error result
            in_range_index=find(temp_fitting_result(:,num_fitted_param+num_SDS_cw+num_SDS_tr+3)<=temp_fitting_result(1,num_fitted_param+num_SDS_cw+num_SDS_tr+3)+error_range_threshold);
            assert(length(in_range_index)==max(in_range_index));
            if length(in_range_index)>max_choose_number
                in_range_index=in_range_index(1:max_choose_number);
            end
            
            choosed_fitted_result=temp_fitting_result(in_range_index,:);
            
            for i=1:length(in_range_index)
                temp_fitted_OP=load(fullfile(temp_target_folder,['fitting_' num2str(orig_init_index(i))],'fitted_mu.txt'));
                temp_fitted_OP=temp_fitted_OP(:,OP_column_arr);
%                 OP_error_arr(:,:,OP_error_counter)=abs(temp_fitted_OP./OP_answer_arr(:,:,target_i)-1);
                OP_error_arr(:,:,OP_error_counter)=temp_fitted_OP./OP_answer_arr(:,:,target_i)-1;
                OP_error_counter=OP_error_counter+1;
            end
        end
    end
end

%% calculate
mean_OP_error_arr=mean(OP_error_arr,3);
std_OP_error_arr=std(OP_error_arr,[],3);
total_mean_OP_error_arr=mean(mean_OP_error_arr,1);
total_std_error_arr=mean(std_OP_error_arr,1);

%% calculate error pdf for each OP
n_hist_interval=100;
OP_pdf_arr={};
OP_pdf_binCenter={};
for i=1:3
    temp_OP_error=reshape(OP_error_arr(:,i,:),1,[]);
    [OP_pdf_arr{i},OP_pdf_binCenter{i}]=hist(temp_OP_error,n_hist_interval);
    OP_pdf_arr{i}=OP_pdf_arr{i}/sum(OP_pdf_arr{i});
end

%% calculate confidence interval for each OP
confidence_arr=[68 95 99.7];  % 68 95 99.7
OP_error_CI={};
for i=1:6
    temp_OP_error=sort(reshape(OP_error_arr(:,i,:),1,[]));
    for ci=1:length(confidence_arr)
        tt=(1-confidence_arr(ci)/100)/2;
        tt=[tt 1-tt];
        tt=round(tt*length(temp_OP_error));
        OP_error_CI{i}(ci,1:2)=temp_OP_error(tt);
    end
end

%% save
if do_use_add_error
    save_file_name='OP_error_arr_Error.mat';
else
    save_file_name='OP_error_arr_noError.mat';
end
save(fullfile(input_dir,'arrangement',fitting_dir,save_file_name),'OP_error_arr','mean_OP_error_arr','std_OP_error_arr','total_mean_OP_error_arr','total_std_error_arr','OP_pdf_arr','OP_pdf_binCenter','confidence_arr','OP_error_CI');

%% plot error distribution
op_name_arr={'\mu_{a,scalp}','\mu_{s,scalp}','\mu_{a,skull}','\mu_{s,skull}','\mu_{a,GM}','\mu_{s,GM}'};
ci_line_style={'-.','--',':'};
tile_index_arr=[1 4 2 5 3 6];
fig=figure('Units','pixels','position',[0 0 1600 800]);
ti=tiledlayout('flow','TileSpacing','compact','Padding','none');
for i=1:6
    nexttile(tile_index_arr(i));
    hist(reshape(OP_error_arr(:,i,:),1,[])*100,40);
    for ci=1:length(confidence_arr)
        for j=1:2
            xline(OP_error_CI{i}(ci,j)*100,ci_line_style{ci},num2str(OP_error_CI{i}(ci,j)*100,'%.2f%%'))
        end
    end
    xlabel([op_name_arr{i} ' error(%)']);
    ylabel('count')
    grid on;
    set(gca,'fontsize',fontSize, 'FontName', 'Times New Roman');
end
if do_use_add_error
    print(fullfile(input_dir,'arrangement',fitting_dir,'OP_error_arr_Error_hist.png'),'-dpng','-r200');
else
    print(fullfile(input_dir,'arrangement',fitting_dir,'OP_error_arr_noError_hist.png'),'-dpng','-r200');
end

% close all;

%% plot
% % op_name_arr={'\mu_{a1}','\mu_{s1}','\mu_{a2}','\mu_{s2}','\mu_{a4}','\mu_{s4}'};
% op_name_arr={'\mu_{a,scalp}','\mu_{s,scalp}','\mu_{a,skull}','\mu_{s,skull}','\mu_{a,GM}','\mu_{s,GM}'};
% tile_index_arr=[1 4 2 5 3 6];
% fig=figure('Units','pixels','position',[0 0 1600 800]);
% % set(fig,'visible','off');
% ti=tiledlayout('flow','TileSpacing','compact','Padding','none');
% for i=1:3
%     nexttile(tile_index_arr(i));
%     shadedErrorBar(fitting_wl,mean_OP_error_arr(:,i)*100,std_OP_error_arr(:,i)*100,'lineprops','-b','patchSaturation',0.33);
%     xlabel('wavelength(nm)');
%     ylabel([op_name_arr{i} ' error(%)']);
%     grid on;
%     set(gca,'fontsize',fontSize, 'FontName', 'Times New Roman');
% end
% if do_use_add_error
%     print(fullfile(input_dir,'arrangement',fitting_dir,'OP_error_arr_Error.png'),'-dpng','-r200');
% else
%     print(fullfile(input_dir,'arrangement',fitting_dir,'OP_error_arr_noError.png'),'-dpng','-r200');
% end
% 
% close all;

disp('Done!');