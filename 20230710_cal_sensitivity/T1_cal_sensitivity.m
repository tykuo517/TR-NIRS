%{
Set the change rate for target layer to test the sensitivity of different SDS and gate

Ting-Yi Kuo
Last update: 2023/07/11
%}

clear; close all;
global tr_net tr_param_range;

baseline=[0.275 0.2 0.3 160 125 150]; % [0.2 0.15 0.5 150 180 150]
changerate_to_exam=[-20 -10 10 20];
subject_name_arr={'KB'}; % ,'ZJ','WH'
interp_folder={'KB_test2_2023-07-17-13-27-58'}; %'KB_test1_2023-07-10-18-40-48','BY_test1_2023-07-10-18-43-48'
model_dir='model_arrange';
title_arr={'\mu_{a,scalp}','\mu_{a,skull}','\mu_{a,GM}','\mu_{s,scalp}','\mu_{s,skull}','\mu_{s,GM}'};
ylim_arr={[-0.2 0.4],[-1 1.5],[-0.5 4],[-0.5 0.5],[-0.05 1],[-0.05 0.05]};
num_SDS=5;
num_gate=10;

%% make the mu table
mu_param_arr=[];
testing_index=1;

%% for baseline
temp_mu_arr=baseline;
temp_param=zeros(1,6);
temp_param(1,:)=temp_mu_arr;
mu_param_arr(testing_index,:)=temp_param;

testing_index=testing_index+1;

%% for changing mu
for l=1:size(baseline,2) % the layer to change
    for del_mus=1:length(changerate_to_exam)
        temp_mu_arr=baseline;
        temp_mu_arr(l)=temp_mu_arr(l)*(1+changerate_to_exam(del_mus)/100);
        temp_param=zeros(1,6);
        temp_param(1,:)=temp_mu_arr;
        mu_param_arr(testing_index,:)=temp_param;

        testing_index=testing_index+1;
    end
end

for i=1:size(mu_param_arr,1)
    mu_param_arr_save(i,:)=[mu_param_arr(i,1) mu_param_arr(i,2) 0.042 mu_param_arr(i,3) mu_param_arr(i,4) mu_param_arr(i,5) 23 mu_param_arr(i,6)];
end
save('OP_sim_sen.txt','mu_param_arr','-ascii','-tabs');


%% Calculate senstivity per subject
dtof=zeros(num_gate,num_SDS,size(mu_param_arr,1));
temp_dtof=zeros(num_gate,num_SDS);

use_ANN=1; % ANN=1, interpolation=0

for sbj=1:length(subject_name_arr)
    %% Get dtof
    if use_ANN==1
        load(fullfile(model_dir,[subject_name_arr{sbj} '_tr_model.mat'])); % net, param_range
        for i=1:size(mu_param_arr,1)
            temp_dtof_=fun_ANN_forward(mu_param_arr(i,:),1);
            for s=1:num_SDS
                temp_dtof(:,s)=temp_dtof_((s-1)*num_gate+1:s*num_gate)';
            end
            dtof(:,:,i,sbj)=temp_dtof;
        end
    else 
        load(['/home/md703/Documents/Ty/20230710_MCX_lookup_table_dataGen_testing_SDS/' interp_folder{sbj} '/all_param_arr.mat']);
        for i=1:size(all_param_arr,1)
            temp_dtof_=all_param_arr(i,9:end);
            for s=1:num_SDS
                temp_dtof(:,s)=temp_dtof_((s-1)*num_gate+1:(s-1)*num_gate+num_gate)';
            end
            dtof(:,:,i,sbj)=temp_dtof;
        end
    end
    
    %% Calculate relative change

    relative_change=zeros(num_gate,num_SDS,size(mu_param_arr,1)-1);

    for i=2:size(mu_param_arr,1)
        relative_change(:,:,i-1)=dtof(:,:,i,sbj)./dtof(:,:,1,sbj)-1;
    end

    % compare SDS
    % f=figure('Position',[0 0 1800 600]);
    % ti=tiledlayout(3,4);%(3,4)
    % for l=1:size(mu_param_arr,1)-1
    %     nexttile;
    %     for s=1:4
    %         plot(1:1:num_gate,relative_change(:,s,l)')
    %         hold on;
    %     end
    % end
    % legend('SDS1','SDS2','SDS3','SDS4','SDS5','SDS6');

    % compare change
    f=figure('Position',[0 0 1920 1080]);
    ti=tiledlayout(num_SDS,6);%(3,4)
    for s=1:num_SDS
        for l=1:size(baseline,2)
            nexttile;
            ind=4*(l-1)+1;
            max_value=max(relative_change(:,:,[ind:ind+3]),[],'all');
            min_value=min(relative_change(:,:,[ind:ind+3]),[],'all');
            for ch=1:4
                plot(1:1:num_gate,relative_change(:,s,ind));
                ind=ind+1;
                hold on;
            end
            xlabel('Time gate');
            ylabel('\DeltaR/R_{b}');
            xticks(1:10);
            xlim([1 10]);
            ylim([min_value max_value]);
            title(['SDS' num2str(s) ',' title_arr{l}]);
    %         ylim(ylim_arr{l});
        end
    end
    leg = legend('-20%','-10%','10%','20%','Orientation','horizontal');
    leg.Layout.Tile = 'south';
    title(ti,subject_name_arr{sbj});
    
    mkdir('results');
    if use_ANN==1
        print(fullfile('results',['relative_change_' subject_name_arr{sbj} '.png']),'-dpng','-r200');
    else
        print(fullfile('results',['relative_change_interp_' subject_name_arr{sbj} '.png']),'-dpng','-r200');
    end


    %% Calculate sensitivity

    for i=1:length(baseline)
        index=1+length(changerate_to_exam)*(i-1);
        for j=1:length(changerate_to_exam)
            sens_temp(:,:,j)=relative_change(:,:,index)./(changerate_to_exam(j)/100);
            index=index+1;
        end
        sens(:,:,i,sbj)=sum(sens_temp,3)./length(changerate_to_exam);
    end


    % plot
    fig=figure('Units','pixels','position',[0 0 1920 500]);%,'position',[0 0 1600 600]
    ti=tiledlayout(1,6);

    for i=1:size(sens,3)
        nexttile;
        for s=1:num_SDS
            plot(1:1:num_gate,sens(:,s,i,sbj));
            hold on;
        end
        xlabel('time gate');
        ylabel('sensitivity');
        xticks(1:10);
        xlim([1 10]);
        title(title_arr(i));
    end
    leg = legend('SDS1','SDS2','SDS3','SDS4','SDS5','Orientation','horizontal');
    leg.Layout.Tile = 'south';

    if use_ANN==1
        print(fullfile('results',['sensitivity_ANN_' subject_name_arr{sbj} '.png']),'-dpng','-r200');
    else
        print(fullfile('results',['sensitivity_interp_' subject_name_arr{sbj} '.png']),'-dpng','-r200');
    end
    
end

%% Calculate mean sensitivity for all subject

sens_mean=mean(sens,4);
sens_std=std(sens,[],4);

% plot
fig=figure('Units','pixels','position',[0 0 1920 500]);%,'position',[0 0 1600 600]
ti=tiledlayout(1,6);

colormap_arr=jet(num_SDS);

for i=1:size(sens,3)
    nexttile;
    for s=1:num_SDS
%         shadedErrorBar(1:1:num_gate,squeeze(sens_mean(:,s,i)),squeeze(sens_std(:,s,i)),'lineprops',{'Color',colormap_arr(s,:),'LineWidth',2},'patchSaturation',0.1);  %'-b'
        plot(1:1:num_gate,squeeze(sens_mean(:,s,i)));
        hold on;
    end
    xlabel('time gate');
    ylabel('sensitivity');
    xticks(1:10);
    xlim([1 10]);
    title(title_arr(i));
end
leg = legend('SDS1','SDS2','SDS3','SDS4','SDS5','Orientation','horizontal');
leg.Layout.Tile = 'south';


if use_ANN==1
    print(fullfile('results',['sensitivity_all.png']),'-dpng','-r200');
else
    print(fullfile('results',['sensitivity_interp_all.png']),'-dpng','-r200');
end


%% Calculate sensitivity ratio for all subject

sens_mean_abs=mean(abs(sens),4);

for i=1:length(baseline)
    sen_ratio_arr(:,:,i)=sens_mean_abs(:,:,i)./sum(sens_mean_abs,3);
end

fig=figure('Units','pixels','position',[0 0 1920 500]); 
ti=tiledlayout(1,6);

for i=1:size(sens_mean_abs,3)
    nexttile;
    for s=1:num_SDS
%         [xx,yy]=polyxpoly(1:1:num_gate,sen_ratio_arr(:,s,i),1:1:num_gate,ones(1,num_gate)*0.05);
        hold on;
        plot(1:1:num_gate,sen_ratio_arr(:,s,i));
%         scatter(xx,yy,'LineWidth',2);
    end
    plot(1:1:num_gate,ones(1,num_gate)*0.05,'LineWidth',2);
    xlabel('time gate');
    ylabel('sensitivity ratio');
    xticks(1:10);
    xlim([1 10]);
    grid on;
    box on;
    title(title_arr(i));
end

leg = legend('SDS1','SDS2','SDS3','SDS4','SDS5','Orientation','horizontal');
leg.Layout.Tile = 'south';

if use_ANN==1
    print(fullfile('results',['sensitivity_ratio.png']),'-dpng','-r200');
else
    print(fullfile('results',['sensitivity_ratio_interp.png']),'-dpng','-r200');
end


    