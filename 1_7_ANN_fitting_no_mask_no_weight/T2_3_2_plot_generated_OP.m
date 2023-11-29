
clc;clear;close all; clearvars -global;

global lambda net param_range 

%% param
error_range_threshold=0.005; % plot the fitting in this distance to the best result
max_choose_number=1; % the max number of fitting

input_dir='test_fitting_2023-10-31-21-57-22'; % please move the fitting folders into this folder first.
subject_name_arr={'KB'};%'KB','WH','WW'
num_anser_to_generate=15; % number of target spec (true answer)

% CW setting
num_SDS_cw=6; % how many SDS are in the target spectrum
SDS_choosed_cw=[1 2 3 4 5 6]; % the SDS chosen to fitted in the target spectrum
SDS_sim_correspond_cw=[1 2 3 4 5 6]; % the SDS index in the simulated spectrum corresponding to each SDS in the target spectrum, SDS_sim_correspond(x)=y means 'measure x = sim y'

% TR setting
num_SDS_tr=5;
SDS_choosed_tr=[3 4 5]; % the SDS chosen to fitted in the target spectrum
SDS_sim_correspond_tr=[1 2 3 4 5]; % the SDS index in the simulated spectrum corresponding to each SDS in the target spectrum, SDS_sim_correspond(x)=y means 'measure x = sim y'

num_gate=10;
gate_choosed=[1 2 3 4 5];
gate_sim_correspond=[1 2 3 4 5 6 7 8 9 10];

num_fitted_param=13; % the number of the fitted parameters

fitting_wl_file='fitting_wl.txt'; % the file in the epsilon folder containing the fitting wavelength
fitting_wl=load(fullfile('epsilon',fitting_wl_file));
lambda=fitting_wl;
fitting_wl_tr=810;

model_dir='model_arrange'; % the folder containing the arranged ANN file


%% Load OP answer

OP_column_arr=[1 2 3 4 7 8]; % the column in OP array of each OP

OP_answer_arr=[];
for i=1:num_anser_to_generate
    temp_OP_arr=load(fullfile(input_dir,'answers',['OP_ans_' num2str(i) '.txt']));
    OP_answer_arr(:,:,i)=interp1(temp_OP_arr(:,1),temp_OP_arr(:,1+OP_column_arr),lambda);
end

param_ans_arr=load(fullfile(input_dir,'answers',['param_answer_arr.txt']));

op_name_arr={'\mu_{a,scalp}','\mu_{s,scalp}','\mu_{a,skull}','\mu_{s,skull}','\mu_{a,GM}','\mu_{s,GM}'};
tile_index_arr=[1 4 2 5 3 6];
fig=figure('Units','pixels','position',[0 0 1600 800]);
ti=tiledlayout('flow','TileSpacing','compact','Padding','none');

colormap_arr=repmat([1 0 0],15,1);
colormap_arr(12,:)=[0 0 1];
% colormap_arr(6,:)=[0 0 1];
% colormap_arr(10,:)=[0 0 1];
% colormap_arr(13,:)=[0 0 1];
% colormap_arr(14,:)=[0 0 1];
% colormap_arr(15,:)=[0 0 1];

for i=1:size(OP_answer_arr,2)
    nexttile(tile_index_arr(i));
    for j=1:size(OP_answer_arr,3)
        plot(lambda,squeeze(OP_answer_arr(:,i,j)),'Color',colormap_arr(j,:));
        hold on
    end
    xlabel('wavelength');
    ylabel('OP');
    title(op_name_arr(i));
end
lgd=legend({'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15'});
lgd.Layout.Tile='south';
print(fullfile(input_dir,'answers','generated_OP.png'),'-dpng','-r200');

% plot param
fig=figure('Units','pixels','position',[0 0 1600 800]);
param_name_arr={'A1','K1','A2','K2','A4','K4',' hc_1','sto2_1','hc_2','sto2_2','hc_4','sto2_4','mel'};
% colormap_arr=jet(size(param_ans_arr,1));

colormap_arr=repmat([1 0 0],15,1);
colormap_arr(12,:)=[0 0 1];
% colormap_arr(6,:)=[0 0 1];
% colormap_arr(10,:)=[0 0 1];
% colormap_arr(13,:)=[0 0 1];
% colormap_arr(14,:)=[0 0 1];
% colormap_arr(15,:)=[0 0 1];

for i=1:size(param_ans_arr,1)
%     for j=1:size(param_ans_arr,2)
        plot(1:1:13,param_ans_arr(i,:),'Color',colormap_arr(i,:));
        xticks(1:1:13)
        xticklabels(param_name_arr)
        hold on
%     end
end

legend({'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15'});

