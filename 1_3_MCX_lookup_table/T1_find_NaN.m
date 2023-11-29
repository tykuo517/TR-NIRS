%{
Find which gate or SDS receives no photons 

Ting-Yi Kuo
Last update: 2023/10/18
Version: 4.41
%}

clear;close all;

num_SDS=5;
num_gate=10;

SDS_to_examine=[4 5];
gate_to_examine=[1];

record_table=cell(num_gate,num_SDS);

sbj_arr = {'KB'};
for sbj = 1
    mus_table = load(fullfile(sbj_arr{sbj},'mus_table.txt'));

    for sim = 1:size(mus_table,1)
        
        % receive no photon (use CV to examine)
%         cv_value=load(fullfile(sbj_arr{sbj},['sim_' num2str(sim)],'SDS_CV_arr_1.txt'));
%         
%         for s=1:num_SDS
%             for g=1:num_gate
%                 if isnan(cv_value(g,s)) 
%                     record_table{g,s}=[record_table{g,s}; mus_table(sim,:)];
%                     record_table{g,s}=unique(record_table{g,s}, 'rows');
%                 end
%             end
%         end
        
        % receive less than 5 photons (use pathlength to examine)
        load(fullfile(sbj_arr{sbj},['sim_' num2str(sim)],'PL_1.mat'));
        
        for s=1:num_SDS
            for g=1:num_gate
                if size(SDS_detpt_arr{g,s},1)<=5
                    record_table{g,s}=[record_table{g,s}; mus_table(sim,:)];
                    record_table{g,s}=unique(record_table{g,s}, 'rows');
                end
            end
        end
        
        save(fullfile('results',['no_photons_SDS' sprintf('%d', SDS_to_examine) '_gate' num2str(gate_to_examine) '_' sbj_arr{sbj} '.mat']),'record_table');

    end    
    fprintf('Done!\n');
    
end