%{
Find 

Ting-Yi Kuo
Last update: 2023/10/18
Version: 4.41
%}

clear;close all;

num_SDS=5;
num_gate=10;

sbj_arr={'KB'};
for sbj=1
    mus_table=load(fullfile(sbj_arr{sbj},'mus_table.txt'));
    cv_arr=zeros(num_gate,num_SDS,length(mus_table));

    for sim=1:size(mus_table,1)
        
        cv_arr(:,:,sim)=load(fullfile(sbj_arr{sbj},['sim_' num2str(sim)],'SDS_CV_arr_1.txt'));
              
%         to_save=[target_mua(1:4) mus_table(sim,1:4) dtof];
%         save(fullfile(sbj_arr{sbj},['DTOF_' num2str(step) '.mat']),'to_save');
% 
%         fprintf(['Finish sim ' num2str(sim) '/' num2str(size(mus_table,1)) '\n']);
    end
    
    %% Observe variation of cv caused by changing mus of scalp, skull, GM
    figure('Units','pixels','position',[0 0 1920 1080]);
    ti=tiledlayout(3,6);
    colormap('jet');
    for s=3:5 %1:num_SDS
        for g=1:6 %6:num_gate
            nexttile;
            scatter3(mus_table(:,1), mus_table(:,2), mus_table(:,4), 50, squeeze(cv_arr(g,s,:))/sqrt(10), 'filled');
            colorbar;
            title(['SDS' num2str(s) ' Gate' num2str(g)]);
            xlabel('scalp');
            ylabel('skull');
            zlabel('GM');
        end
    end
    saveas(gcf,fullfile('results',['cv_variation_' sbj_arr{sbj} '.fig']));
    
    %% Choose the cv to use for generating target with noise (in 1_7_ANN_fitting)
    choosed_index=336; % [150 125 23 150]: middle of all mus combinations
    cv_value=cv_arr(:,:,choosed_index)/sqrt(10);
    save(fullfile('results','cv_1E11.mat'),'cv_value');
        
    fprintf('Done!\n');
    
end
