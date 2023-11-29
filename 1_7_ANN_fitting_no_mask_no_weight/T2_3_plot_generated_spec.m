%{
Plot the generated target spectrum

Benjamin Kao
Last update: 2021/01/17
%}

clc;clear;close all;

global lambda Lbound Ubound net param_range;

%% param 
subject_name_arr={'KB'};
num_anser_to_generate=15; % number of target spec (true answer)
num_error_to_generate=15; % number of adding noise to the same, the first one will have no error
num_SDS_cw=6;

num_SDS_tr=5;
num_gate=10;
SDS_dist_arr_tr=[1.5 2.2 2.9 3.6 4.3]; % cm
SDS_dist_arr_cw=[0.8 1.5 2.12 3 3.35 4.5 4.74]; % cm; % cm
input_dir='test_fitting_2023-10-31-21-57-22';

%% main
for sbj=1:length(subject_name_arr)
    for target_i=1:num_anser_to_generate
        spec_arr=[];
        for error_i=1:num_error_to_generate
            spec_arr(:,:,error_i)=load(fullfile(input_dir,subject_name_arr{sbj},['target_' num2str(target_i) '_' num2str(error_i)],'target_spec.txt')); % CW
            dtof_arr(:,:,error_i)=load(fullfile(input_dir,subject_name_arr{sbj},['target_' num2str(target_i) '_' num2str(error_i)],'target_dtof.txt')); % TR
        end
        
        % CW
        mean_spec=mean(spec_arr,3);
        std_spec=std(spec_arr,[],3);
        cv_spec=std_spec./mean_spec;
        
        % TR
        mean_dtof=mean(dtof_arr,3);
        std_dtof=std(dtof_arr,[],3);
        cv_dtof=std_dtof./mean_dtof;        
        
        
        % CW
        fig=figure('Position',[0 0 1920 1080]);
        set(fig,'visible','off');
        ti=tiledlayout('flow','TileSpacing','compact','Padding','compact');
        for s=1:num_SDS_cw
            nexttile();
            plot(spec_arr(:,1,1),squeeze(spec_arr(:,s+1,:)),'LineWidth',2);
            title(['SDS = ' num2str(SDS_dist_arr_cw(s)) ' cm, CV=' num2str(cv_spec(1,s+1)*100,'%.2f%%')]);
            yylim=ylim;
            yylim(1)=0;
            ylim(yylim);
            set(gca,'fontsize',17, 'FontName', 'Times New Roman');
            grid on;
        end
        title(ti,[subject_name_arr{sbj} ' target ' num2str(target_i)]);
        print(fullfile(input_dir,subject_name_arr{sbj},['target_spec_' num2str(target_i) '.png']),'-dpng','-r200');

        % TR
        fig=figure('Position',[0 0 1920 1080]);
        set(fig,'visible','off');
        ti=tiledlayout('flow','TileSpacing','compact','Padding','compact');
        for s=1:num_SDS_tr
            nexttile();
            semilogy(1:1:10,squeeze(dtof_arr(:,s,:)),'LineWidth',1);
%             title(['SDS = ' num2str(SDS_dist_arr(s)) ' cm, CV=' num2str(cv_spec(1,s+1)*100,'%.2f%%')]);
            title(['SDS = ' num2str(SDS_dist_arr_tr(s)) ' cm']);
            yylim=ylim;
            yylim(1)=0;
            ylim(yylim);
            ylabel('Reflectance');
            yyaxis right
            plot(1:1:10,cv_dtof(:,s)*100);
            ylabel('CV (%)');
            set(gca,'fontsize',10, 'FontName', 'Times New Roman');
            grid on;
        end
        title(ti,[subject_name_arr{sbj} ' target ' num2str(target_i)],'fontsize',10);
        print(fullfile(input_dir,subject_name_arr{sbj},['target_dtof_' num2str(target_i) '.png']),'-dpng','-r200');
        close all;
    end
end

disp('Done!');