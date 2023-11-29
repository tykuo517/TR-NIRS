%{
Split the pathlength of simulations into time gates, calculate DTOF of each simulation, 
and store in each simulation folder

Ting-Yi Kuo
Last update: 2023/3/29
Version: 4.41
%}

% clear;close all;


num_SDS=7;

sbj_arr = {'KB'};
for sbj = 1
    mus_table = load(fullfile(sbj_arr{sbj},'mus_table.txt'));
%     step=1;

    for sim = 1:size(mus_table,1)
        %% Calculate reflectance
        load(fullfile(sbj_arr{sbj},['sim_' num2str(sim)],'PL_1.mat'));
        
        target_mua=[0.6 0.45 0.1 0.5];
        target_mua(:,5)=target_mua(:,4)*0.5;
        target_mua(:,6)=0;
        
        reflectance=zeros(1,num_SDS);
        index=1;
        for s=1:num_SDS
            if size(SDS_detpt_arr{s},1)>0
                reflectance(index)=1/each_photon_weight_arr(s)*sum(exp(-1*sum(double(SDS_detpt_arr{s}).*target_mua,2)),1);%*(true_r/sim_set.detector_r).^2;
            else
                reflectance(index)=0;
            end
            index=index+1;
        end
        
        
        to_save=[target_mua(1:4) mus_table(sim,1:4) reflectance];
        save(fullfile(sbj_arr{sbj},['reflectance_' num2str(sim) '.mat']),'to_save');
%         step=step+1;

        fprintf(['Finish sim ' num2str(sim) '/' num2str(size(mus_table,1)) '\n']);
    end
        
    fprintf('Done!\n');
    
end


    