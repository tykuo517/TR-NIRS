%{
Compare the fitted spec simulated using ANN and MC, also compare with the target spec

Benjamin Kao
Last update: 2021/04/05
%}

clc;clear;close all;

%% param
ANN_dir='fitted_result';
MC_dir='~/Documents/BenjaminKao/20210220_simulate_fitted_OP/pathlength/multiSol/sim_5E9_n1457_diffNA_Andy/cal_PL_allSDS';

output_dir='compare_ANN_MC';

target_dir='input_target_2'; % the folder containing the target spectrum
target_name_arr={'tc_mean','ww_mean','wh_mean','yf_mean','kb_mean'}; % the name of the target spectrum
target_MC_sol_index=[2 1 1 1 1]; % because the MC dir has multisol, so should choose which result to use
num_SDS=6;
sim_detector_r=[1 2 2 2 2 2]; % mm
true_detector_r=[0.2 0.2 0.2 0.2 0.2 0.2]; % mm
SDS_distance_arr=[0.8 1.5 2.12 3 3.35 4.5]; % cm

%% init
if exist(fullfile(ANN_dir,output_dir),'dir')==0
    mkdir(fullfile(ANN_dir,output_dir));
end

SDS_ratio_arr=(sim_detector_r./true_detector_r).^2; % the ratio between simulated r and true r

%% main
sbj_SDS_error_arr=[];

for sbj=1:length(target_name_arr)
    target_spec=load(fullfile(target_dir,[target_name_arr{sbj} '.txt']));
    ANN_spec=load(fullfile(ANN_dir,[target_name_arr{sbj} '_fitted_spec_1.txt']));
    MC_spec=load(fullfile(MC_dir,[target_name_arr{sbj} '_' num2str(target_MC_sol_index(sbj)) '_reflectance.txt']));
    MC_spec(:,2:1+num_SDS)=MC_spec(:,2:1+num_SDS)./SDS_ratio_arr;
    interp_wl=MC_spec(find(MC_spec(:,1)>=min(ANN_spec(:,1)) & MC_spec(:,1)<=max(ANN_spec(:,1))),1);
    ANN_spec_interp=interp1(ANN_spec(:,1),ANN_spec(:,2:end),interp_wl);
    MC_spec_interp=interp1(MC_spec(:,1),MC_spec(:,2:end),interp_wl);
    target_spec_interp=interp1(target_spec(:,1),target_spec(:,2:end),interp_wl);
    target_ANN_error=sqrt(mean((ANN_spec_interp(:,1:num_SDS)./target_spec_interp-1).^2,1));
    target_MC_error=sqrt(mean((MC_spec_interp(:,1:num_SDS)./target_spec_interp-1).^2,1));
    SDS_error_arr=sqrt(mean((MC_spec_interp./ANN_spec_interp-1).^2,1));
    sbj_SDS_error_arr(sbj,:)=SDS_error_arr;
    
%     fig=figure('Units','inches','position',[0 0 7.165 4.5]);
    fig=figure('Units','inches','position',[0 0 14.33 9]);
    ti=tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
    for s=1:num_SDS
        nexttile();
        plot(target_spec(:,1),target_spec(:,1+s),ANN_spec(:,1),ANN_spec(:,1+s),MC_spec(:,1),MC_spec(:,1+s),'LineWidth',2);
        legend({'target',['ANN, err=' num2str(target_ANN_error(s)*100,'%.2f%%')],['MC, err=' num2str(target_MC_error(s)*100,'%.2f%%')]},'Location','southeast');
        xlabel('wavelength(nm)');
        ylabel('reflectance');
        xlim([min(ANN_spec(:,1)),max(ANN_spec(:,1))]);
        ylim([0 inf]);
        grid on;
%         title(['SDS = ' num2str(SDS_distance_arr(s)) ' cm, error = ' num2str(SDS_error_arr(s)*100,'%.2f%%')],'fontsize',18, 'FontName', 'Times New Roman');
        title(['SDS = ' num2str(SDS_distance_arr(s)) ' cm'],'fontsize',18, 'FontName', 'Times New Roman');
        set(gca,'fontsize',12, 'FontName', 'Times New Roman');
    end
    title(ti,['Subject ' num2str(sbj) ' fitted spectrum'],'fontsize',22, 'FontName', 'Times New Roman')
    print(fullfile(ANN_dir,output_dir,['subject_' num2str(sbj) '.png']),'-dpng','-r300');
    close all;
end

save(fullfile(ANN_dir,output_dir,'subject_SDS_error_arr.txt'),'sbj_SDS_error_arr','-ascii','-tabs');

disp('Done!');