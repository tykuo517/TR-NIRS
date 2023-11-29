%{
Set the change rate for target layer to test the sensitivity of different SDS and gate

Ting-Yi Kuo
Last update: 2023/07/11
%}

clear; %close all;
global tr_net tr_param_range;

gate_checked=[3 5 7];
plot_individual=1;

baseline=[0.15 0.15 0.2 120 75 125]; %  [0.275 0.2 0.3 160 125 150] [0.2 0.3 0.2 150 125 150] today:[0.15 0.25 0.3 150 125 150] [0.275 0.2 0.25 160 125 150] [0.15 0.15 0.2 120 75 125] [0.35 0.25 0.3 200 175 175]
changerate_to_exam=[-20 -10 10 20];
subject_name_arr={'WH'}; % ,'ZJ','WH'
interp_folder={'KB_test2_2023-07-17-13-27-58'}; %'KB_test1_2023-07-10-18-40-48','BY_test1_2023-07-10-18-43-48'
model_dir='model_arrange';
title_arr={'\mu_{a,scalp}','\mu_{a,skull}','\mu_{a,GM}','\mu_{s,scalp}','\mu_{s,skull}','\mu_{s,GM}'};
SDS_arr=[1.5 2.2 2.9 3.6 4.3];
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
dtof=zeros(length(gate_checked),num_SDS,size(mu_param_arr,1));
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
            dtof(:,:,i,sbj)=temp_dtof(gate_checked,:);
        end
    else 
        load(['/home/md703/Documents/Ty/20230710_MCX_lookup_table_dataGen_testing_SDS/' interp_folder{sbj} '/all_param_arr.mat']);
        for i=1:size(all_param_arr,1)
            temp_dtof_=all_param_arr(i,9:end);
            for s=1:num_SDS
                temp_dtof(:,s)=temp_dtof_((s-1)*num_gate+1:(s-1)*num_gate+num_gate)';
            end
            dtof(:,:,i,sbj)=temp_dtof(gate_checked,:);
        end
    end
    
%     dtof=dtof(gate_checked,:,:,:);
    

    %% Plot reflectance change
    if plot_individual
%         for g=1:length(gate_checked)
%             f=figure('Position',[0 0 1920 1080]);
%             set(f,'visible','off');
%             ti=tiledlayout(6,num_SDS);%(3,4)
%             for l=1:size(baseline,2)
%                 for s=1:num_SDS
%                     nexttile;
%                     ind=4*(l-1)+2;
% 
%                     plot(changerate_to_exam,squeeze(dtof(g,s,ind:ind+3,sbj)),'-o');
% 
%                     xlabel(title_arr(l));
%                     ylabel('reflectance');
%                     title(['SDS' num2str(s)]);
%                 end
%             end
%             title(ti,['Gate ' num2str(gate_checked(g)) ' ' subject_name_arr{sbj}]);
% 
%             if use_ANN==1
%                 print(fullfile('results',['gate' num2str(gate_checked(g)) '_reflectance_' subject_name_arr{sbj} '.png']),'-dpng','-r200');
%             else
%                 print(fullfile('results',['relative_change_interp_' subject_name_arr{sbj} '.png']),'-dpng','-r200');
%             end
%         end
    end
    

    %% Calculate sensitivity
    
    % calculate relative change
    relative_change=zeros(length(gate_checked),num_SDS,size(mu_param_arr,1)-1);

    for i=2:size(mu_param_arr,1)
        relative_change(:,:,i-1)=dtof(:,:,i,sbj)./dtof(:,:,1,sbj)-1;
    end
    
    % calculate sensitivity
    for i=1:length(baseline)
        index=1+length(changerate_to_exam)*(i-1);
        for j=1:length(changerate_to_exam)
            sens_temp(:,:,j)=relative_change(:,:,index)./(changerate_to_exam(j)/100);
            index=index+1;
        end
        sens(:,:,i,sbj)=sum(sens_temp,3)./length(changerate_to_exam);
    end
    
    % plot relative change
    if plot_individual
%         for g=1:length(gate_checked)
%             f=figure('Position',[0 0 1920 1080]);
%             set(f,'visible','off');
%             ti=tiledlayout(6,num_SDS);%(3,4)
%             for l=1:size(baseline,2)
%                 for s=1:num_SDS
%                     nexttile;
%                     ind=4*(l-1)+1;
%                     plot(changerate_to_exam,squeeze(relative_change(g,s,ind:ind+3)),'-o');
% 
%                     xlabel(title_arr(l));
%                     ylabel('relative change');
%                     title(['SDS' num2str(s)]);
%                 end
%             end
%             title(ti,['Gate ' num2str(gate_checked(g)) ' ' subject_name_arr{sbj}]);
% 
%             if use_ANN==1
%                 print(fullfile('results',['gate' num2str(gate_checked(g)) '_relative_change_' subject_name_arr{sbj} '.png']),'-dpng','-r200');
%             else
%                 print(fullfile('results',['relative_change_interp_' subject_name_arr{sbj} '.png']),'-dpng','-r200');
%             end
%         end


        % plot sensitivity per subject
        fig=figure('Units','pixels','position',[0 0 1280 720]);%,'position',[0 0 1600 600]
        ti=tiledlayout(2,3);
        for i=1:size(sens,3)
            nexttile;
            for g=1:length(gate_checked)
                plot(SDS_arr,squeeze(sens(g,:,i,sbj)),'Linewidth',2);
                hold on
            end
            xlabel('SDS (cm)');
            ylabel('sensitivity');
            title(title_arr(i));
        end
        title(ti,[subject_name_arr{sbj}]);
        lgd=legend('Gate 3','Gate 5','Gate 7','Orientation','horizontal');
        lgd.Layout.Tile='south';
    end
    
end

%% Calculate mean sensitivity for all subject

sens_mean=mean(sens,4);
sens_std=std(sens,[],4);

% % plot separately
% colormap_arr=jet(num_SDS);
% 
% for g=1:length(gate_checked)
%     fig=figure('Units','pixels','position',[0 0 1920 1080]);%,'position',[0 0 1600 600]
%     ti=tiledlayout(2,3);
%     for i=1:size(sens,3)
%         nexttile;
%         shadedErrorBar(SDS_arr,squeeze(sens_mean(g,:,i)),squeeze(sens_std(g,:,i)),'lineprops',{'LineWidth',2},'patchSaturation',0.1);  %'-b'
%     %         plot(1:1:num_gate,squeeze(sens_mean(:,s,i)));
%         xlabel('SDS (cm)');
%         ylabel('sensitivity');
%         title(title_arr(i));
%     end
%     title(ti,['Gate ' num2str(gate_checked(g))]);
% end

% plot together
colormap_arr=jet(num_SDS);

fig=figure('Units','pixels','position',[0 0 1920 1080]);%,'position',[0 0 1600 600]
ti=tiledlayout(2,3);

for i=1:size(sens,3)
    nexttile;
    for g=1:length(gate_checked)
        shadedErrorBar(SDS_arr,squeeze(sens_mean(g,:,i)),squeeze(sens_std(g,:,i)),'lineprops',{'LineWidth',2},'patchSaturation',0.1);  %'-b'
    %         plot(1:1:num_gate,squeeze(sens_mean(:,s,i)));
        hold on
        grid on
    end
    xlabel('SDS (cm)');
    ylabel('sensitivity (%reflectance change/%OP change)');
    title(title_arr(i));
end
lgd=legend('Gate 3','Gate 5','Gate 7','Orientation','horizontal');
lgd.Layout.Tile = 'south';

if use_ANN==1
    print(fullfile('results','gate_sensitivity.png'),'-dpng','-r200');
else
    print(fullfile('results',['relative_change_interp_' subject_name_arr{sbj} '.png']),'-dpng','-r200');
end


%% Calculate sensitivity ratio for all subject

sens_mean_abs=mean(abs(sens),4);

for i=1:length(baseline)
    sen_ratio_arr(:,:,i)=sens_mean_abs(:,:,i)./sum(sens_mean_abs,3);
end

fig=figure('Units','pixels','position',[0 0 1920 1080]); 
ti=tiledlayout(2,3);

colormap_arr=winter(length(gate_checked));

for i=1:size(sens_mean_abs,3)
    nexttile;
    for g=1:length(gate_checked)

%         [xx,yy]=polyxpoly(1:1:num_gate,sen_ratio_arr(:,s,i),1:1:num_gate,ones(1,num_gate)*0.05);
        hold on;
        if sum(sen_ratio_arr(g,:,i)<0.05)==0
            line_prop='-o';
        elseif sum(sen_ratio_arr(g,:,i)>0.05)==0
            line_prop='-x';
        else
            line_prop='-^';
        end
            
        plot(SDS_arr,sen_ratio_arr(g,:,i),line_prop,'Linewidth',1.5,'Color',colormap_arr(g,:),'HandleVisibility','off');
        plot(SDS_arr,sen_ratio_arr(g,:,i),'Linewidth',1.5,'Color',colormap_arr(g,:));
        
    end
    plot(SDS_arr,ones(1,num_SDS)*0.05,'--','LineWidth',1.5);
    xlabel('SDS (cm)');
    ylabel('sensitivity ratio');
    xlim([SDS_arr(1) SDS_arr(end)]);
    grid on;
    box on;
    title(title_arr(i));
end

leg = legend('Gate 3','Gate 5','Gate 7','Orientation','horizontal');
leg.Layout.Tile = 'south';

if use_ANN==1
    print(fullfile('results',['gate_sensitivity_ratio.png']),'-dpng','-r200');
else
    print(fullfile('results',['gate_sensitivity_ratio_interp.png']),'-dpng','-r200');
end