%{
Plot the fitting result, and compare to the literature OP

Benjamin Kao
Last update: 2021/03/17
%}

clc;clear;close all;

global lambda net param_range 

%% param
input_dir='fitted_result'; % the folder of the fitting result
OP_dir='literature_OPs'; % the folder containing the literature OP files

to_plot_layer=[1 1 0 1]; % if plot scalp, skull, CSF and GM
add_OP_error=1; % if =1, add the OP error to the fitted result
use_CI_as_error=1; % if =1, use confidence interval rather than CI to plot
CI_index=2; % if use_CI_as_error, the confidence interval to plot

% which OP for each type to plot
% 700~900
% do_plot_arr={[1 1 0 0 1 1 0],[0 1 1 1 1 1],[1 1],[1 1 1 1 0 1 1 1]
%     [1 1 0 1 1 1],[1 1 1 1],[1 1],[1 1 1 1 1 1 1]};
% 650~1000
do_plot_arr={[1 1 0 0 1 1 0],[0 1 1 1 1 1],[1 1],[1 1 1 1 0 1 1 1]
    [1 1 0 1 1 1],[1 1 1 1],[1 1],[1 1 1 1 1 1 1]};

target_name_arr={'tc_mean','ww_mean','wh_mean','yf_mean','kb_mean'}; % the name of the target spectrum

%% init
layer_name_arr={'scalp','skull','cerebral spinal fliud','gray matter'};
OP_name={'\mu_a','\mu_s'''};

lit_OP_subdir_name={'mua1_mm','mua2_mm','mua3_mm','mua4_mm'
    'musp1_mm','musp2_mm','musp3_mm','musp4_mm'};

%% main
% load the fitted OP
fitted_OP_arr={};
subject_legend_arr={};
sbj_OP_CV_arr=[];
sbj_OP_CI_arr={};
result_subject_arr=[]; % the subject index of each result
result_multiSolution_arr=[]; % the index of solution of each result
for i=1:length(target_name_arr)
    choosed_rank=load(fullfile(input_dir,[target_name_arr{i} '_choosed_rank.txt']));
    for j=1:length(choosed_rank)
        result_subject_arr(end+1)=i;
        result_multiSolution_arr(end+1)=j;
        subject_legend_arr{end+1}=['subject' num2str(i) '_' num2str(j)];
        fitted_OP_arr{end+1}=load(fullfile(input_dir,[target_name_arr{i} '_fitted_OP_' num2str(j) '.txt']));
        sbj_OP_CV_arr(end+1,:)=load(fullfile(input_dir,[target_name_arr{i} '_OP_CV.txt']));
        temp_OP_CI=load(fullfile(input_dir,[target_name_arr{i} '_error_CI.mat']));
        sbj_OP_CI_arr{end+1}=temp_OP_CI.error_CI;
    end
end

%% plot the OP
multiSolution_style_arr={'-','-.','--'};

OP_index=1;
for layer_i=1:4
    if to_plot_layer(layer_i)
        fig=figure('Units','inches','position',[0 0 7.165 4.5]);
        ti=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
        
        for OP_i=1:2
            OP_input_dir=fullfile(OP_dir,lit_OP_subdir_name{OP_i,layer_i});

            OP_file_arr=jsondecode(fileread(fullfile(OP_input_dir,'OP_file_arr.json')));
            for i=1:length(OP_file_arr.filename_arr)
                fprintf('File %d: %s\n',i,OP_file_arr.filename_arr{i});
            end
            assert(length(do_plot_arr{OP_i,layer_i})==length(OP_file_arr.filename_arr));

            OP_spec_arr={};
            OP_legend_arr={};
            for i=1:length(OP_file_arr.filename_arr)
                if do_plot_arr{OP_i,layer_i}(i)
                    OP_spec_arr{end+1,1}=load(fullfile(OP_input_dir,OP_file_arr.filename_arr{i}));
                    OP_legend_arr{end+1,1}=OP_file_arr.file_legend_arr{i};
                end
            end
            for i=1:length(subject_legend_arr)
                OP_legend_arr{end+1,1}=subject_legend_arr{i};
            end
            
            color_map_arr=lines(length(OP_spec_arr));
            subject_color_map_arr=lines(length(target_name_arr));
            
            nexttile(OP_i);
            hold on;
            % plot literature OP
            for i=1:length(OP_spec_arr)
                plot(OP_spec_arr{i}(:,1),OP_spec_arr{i}(:,2),':','Color',color_map_arr(i,:),'LineWidth',1.5);
            end
            
            if OP_i==1
                scale_multiplyer=1/10; % turn mua from 1/cm into 1/mm
            else
                scale_multiplyer=1/10*(1-0.9); % turn mus 1/cm into musp 1/mm
            end
            
            opii=layer_i*2+OP_i-2;
            if layer_i==4
                opii=opii-2;
            end
            
            % plot the fitted OP
            for sbj=1:length(fitted_OP_arr)
                if add_OP_error==0
                    plot(fitted_OP_arr{sbj}(:,1),fitted_OP_arr{sbj}(:,-1+2*layer_i+OP_i)*scale_multiplyer,multiSolution_style_arr{result_multiSolution_arr(sbj)},'Color',subject_color_map_arr(result_subject_arr(sbj),:),'LineWidth',1.5);
                elseif add_OP_error==1
                    if use_CI_as_error
                        shadedErrorBar(fitted_OP_arr{sbj}(:,1),fitted_OP_arr{sbj}(:,-1+2*layer_i+OP_i)*scale_multiplyer,fitted_OP_arr{sbj}(:,-1+2*layer_i+OP_i)'*scale_multiplyer.*abs(sbj_OP_CI_arr{sbj}{opii}(CI_index,2:-1:1))','lineProps',{multiSolution_style_arr{result_multiSolution_arr(sbj)},'Color',subject_color_map_arr(result_subject_arr(sbj),:),'LineWidth',1.5},'patchSaturation',0.15);
                    else
                        shadedErrorBar(fitted_OP_arr{sbj}(:,1),fitted_OP_arr{sbj}(:,-1+2*layer_i+OP_i)*scale_multiplyer,fitted_OP_arr{sbj}(:,-1+2*layer_i+OP_i)*scale_multiplyer*sbj_OP_CV_arr(sbj,OP_index)*0.01,'lineProps',{multiSolution_style_arr{result_multiSolution_arr(sbj)},'Color',subject_color_map_arr(result_subject_arr(sbj),:),'LineWidth',1.5},'patchSaturation',0.15);                        
                    end
                end
            end
            
            xlim([min(fitted_OP_arr{1,1}(:,1)) max(fitted_OP_arr{1,1}(:,1))]);
            yylim=ylim(); yylim(1)=0; ylim(yylim);
            xlabel('wavelength(nm)');
            ylabel([OP_name{OP_i} ' (1/mm)']);
%             set(gca,'YScale','log');
            grid on;
            lgd=legend(OP_legend_arr,'Location','southoutside');
            lgd.NumColumns = 2;
            set(gca,'fontsize',9, 'FontName', 'Times New Roman');
            OP_index=OP_index+1;
        end
        title(ti,layer_name_arr{layer_i},'fontsize',20, 'FontName', 'Times New Roman');
        
        if add_OP_error==0
            print(fullfile(input_dir,['litOP_L' num2str(layer_i) '.png']),'-dpng','-r600');
        elseif add_OP_error==1
            if use_CI_as_error
                print(fullfile(input_dir,['litOP_L' num2str(layer_i) '_shaded_CI' num2str(CI_index) '.png']),'-dpng','-r600');
            else
                print(fullfile(input_dir,['litOP_L' num2str(layer_i) '_shaded.png']),'-dpng','-r600');
            end
        end
        close all;
    end
end

disp('Done!');