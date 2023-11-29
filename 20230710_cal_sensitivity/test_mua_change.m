%{
Set the change rate for target layer to test the sensitivity of different SDS and gate

Ting-Yi Kuo
Last update: 2023/07/11
%}

clear; close all;
global net param_range;

baseline=[0.2 0.15 245];
changerate_to_exam=[-20 -10 0 10 20];
subject_name_arr={'KB','WH','WW'}; % 'KB','WH','WW'
interp_folder={'KB_test2_2023-07-17-13-27-58'}; %'KB_test1_2023-07-10-18-40-48','BY_test1_2023-07-10-18-43-48'
model_dir='model_arrange';
title_arr={'\mu_{a,skull}','\mu_{a,GM}','\mu_{s,GM}'};
ylim_arr={[-1 1],[-1 5],[-0.5 0.5]};
num_SDS=5;
num_gate=10;

%% make the mu table
mu_param_arr=[];
testing_index=1;

%% for baseline
temp_mu_arr=baseline;
temp_param=zeros(1,3);
temp_param(1,:)=temp_mu_arr;
mu_param_arr(testing_index,:)=temp_param;

testing_index=testing_index+1;

%% for changing mu
for l=2 %:size(baseline,2) % the layer to change
    for del_mus=1:length(changerate_to_exam)
        temp_mu_arr=baseline;
        temp_mu_arr(l)=temp_mu_arr(l)*(1+changerate_to_exam(del_mus)/100);
        temp_param=zeros(1,3);
        temp_param(1,:)=temp_mu_arr;
        mu_param_arr(testing_index,:)=temp_param;

        testing_index=testing_index+1;
    end
end

for i=1:size(mu_param_arr,1)
    mu_param_arr_save(i,:)=[0.46 mu_param_arr(i,1) 0.042 mu_param_arr(i,2) 55 64 23 mu_param_arr(i,3)];
end
save('OP_sim_sen.txt','mu_param_arr_save','-ascii','-tabs');

%% get dtof

dtof=zeros(num_gate,num_SDS,size(mu_param_arr,1));
temp_dtof=zeros(num_gate,num_SDS);

% ANN
use_ANN=1;
for sbj=1:length(subject_name_arr)
    load(fullfile(model_dir,[subject_name_arr{sbj} '_model.mat'])); % net, param_range
    for i=1:size(mu_param_arr,1)
        temp_dtof_=fun_ANN_forward(mu_param_arr(i,:));
        for s=1:num_SDS
            temp_dtof(:,s)=temp_dtof_((s-1)*num_gate+1:(s-1)*num_gate+num_gate)';
        end
        dtof(:,:,i,sbj)=temp_dtof;
    end
end
mean_dtof=mean(dtof,4);

% LUT
% use_ANN=0;
% for sbj=1:length(interp_folder)
%     load(['/home/md703/Documents/Ty/20230710_MCX_lookup_table_dataGen_testing_SDS/' interp_folder{sbj} '/all_param_arr.mat']);
%     for i=1:size(all_param_arr,1)
%         temp_dtof_=all_param_arr(i,9:end);
%         for s=1:num_SDS
%             temp_dtof(:,s)=temp_dtof_((s-1)*num_gate+1:(s-1)*num_gate+num_gate)';
%         end
%         dtof(:,:,i,sbj)=temp_dtof;
%     end
% end
% mean_dtof=mean(dtof,4);


%% Calculate relative change

relative_change=zeros(num_gate,num_SDS,size(mu_param_arr,1)-1);

for i=2:size(mu_param_arr,1)
    relative_change(:,:,i-1)=mean_dtof(:,:,i)./mean_dtof(:,:,1)-1;
end

% compare change
for s=1:num_SDS
    f=figure('Position',[0 0 1920 1080]);
    ti=tiledlayout("flow");
    for g=1:num_gate
        nexttile;
        plot(1:1:length(changerate_to_exam),squeeze(relative_change(g,s,:)));
        xticks(1:9);
        xticklabels({'-20', '-10', '0', '10', '20'});
        xlabel('OP change(%)');
        ylabel('\DeltaR/R_{b}');
        title(['gate' num2str(g)]);
    end
    title(ti,['SDS' num2str(s) ', ' title_arr{l}]);
    ylim(ylim_arr{l});

end

mkdir('results');
if use_ANN==1
    print(fullfile('results','relative_change_SDS_ANN.png'),'-dpng','-r200');
else
    print(fullfile('results','relative_change_SDS_interp.png'),'-dpng','-r200');
end




    