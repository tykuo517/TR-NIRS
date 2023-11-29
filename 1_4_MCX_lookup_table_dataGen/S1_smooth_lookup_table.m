%{
Smooth the mcx lookup table with designated mua parameter set

This code will save files and figures in ./smooth_result/{subject_name}:
'lkt_ref_value_arr_bs.mat'-look up table value before smoothing
'lkt_ref_value_arr_as.mat'-look up table value after smoothing
'mua_param_arr.mat'-mua set used for smoothing 
'in_place_arr'

Ting-Yi Kuo
Last update: 2023/10/16
%}

clc;clear;close all;

%% param
lookup_table_arr='../1_3_MCX_lookup_table'; % the dir containing the unmerged lookup table
subject_name_arr={'WH'}; % the name of the subjects


% mus boundary around 800nm  
% mus_ub=[225 200 23 250]; % 1/cm
% mus_lb=[100 50 23 50]; % 1/cm

% actual simulation of mus boundary
mus_ub=[250 225 275]; % 1/cm, skip CSF
mus_lb=[75 25 25]; % 1/cm, skip CSF

mua_ub=[0.45 0.3 0.042 0.4]; % 1/cm
mua_lb=[0.1 0.1 0.042 0.1]; % 1/cm

num_SDS=5;
num_gate=10;
num_total=num_SDS*num_gate;

num_layer=4; % number of layer to random

max_mua_sameTime=10; % how many mua set to calculate at the same time, use smaller value for smaller memory consumption
linewidth=1.25;

for sbj_i=1:length(subject_name_arr)
    
    subject_name=subject_name_arr{sbj_i};

    %% init
    lkt_dir=fullfile(lookup_table_arr,subject_name); % the dir containing the unmerged lookup table
    output_dir=fullfile('results_smooth',subject_name);
    mkdir(output_dir);

    % load lookup table information
    lkt_sim_set=load(fullfile(lkt_dir,'sim_set.mat'));
    lkt_sim_set=lkt_sim_set.sim_set;
    lkt_layer_mus=load(fullfile(lkt_dir,'layer_mus.mat'));
    lkt_layer_mus=lkt_layer_mus.layer_mus;
    lkt_mus_table=load(fullfile(lkt_dir,'mus_table.txt'));

    % save the bound
    param_ub=[mua_ub mus_ub];
    param_lb=[mua_lb mus_lb];
    to_save=[param_ub; param_lb];
    
    %% Decide your mua set (use lookup table mus as mus set here)
    lkt_layer_mua={[0.1:0.035:0.45],[0.1:0.02:0.3],0.042,[0.1:0.03:0.4]}; % mua for each layer, 1/cm, generate 11x11x1x11=1331
    mua_param_arr=[];
    for i=1:length(lkt_layer_mua{1})
        for j=1:length(lkt_layer_mua{2})
            for k=1:length(lkt_layer_mua{3})
                for l=1:length(lkt_layer_mua{4})
                    mua_param_arr=[mua_param_arr; lkt_layer_mua{1}(i) lkt_layer_mua{2}(j) lkt_layer_mua{3}(k) lkt_layer_mua{4}(l) lkt_layer_mua{4}(l)*0.5];
                end
            end
        end
    end
    
    


    %% find the index in the lookup table, which is a 4-D array, which dimention is the mus for one layer, and each value is the corresponding index in a 1-D array
    in_place_arr=zeros(length(lkt_layer_mus{1}),length(lkt_layer_mus{2}),length(lkt_layer_mus{3}),length(lkt_layer_mus{4}));
    for L1=1:length(lkt_layer_mus{1})
        for L2=1:length(lkt_layer_mus{2})
            for L3=1:length(lkt_layer_mus{3})
                for L4=1:length(lkt_layer_mus{4})
                    lkt_index=find(lkt_mus_table(:,1)==lkt_layer_mus{1}(L1) & lkt_mus_table(:,2)==lkt_layer_mus{2}(L2) & lkt_mus_table(:,3)==lkt_layer_mus{3}(L3) & lkt_mus_table(:,4)==lkt_layer_mus{4}(L4));
                    in_place_arr(L1,L2,L3,L4)=lkt_index;
                end
            end
        end
    end

    in_place_arr_mua=zeros(length(lkt_layer_mua{1}),length(lkt_layer_mua{2}),length(lkt_layer_mua{3}),length(lkt_layer_mua{4}));
    for L1=1:length(lkt_layer_mua{1})
        for L2=1:length(lkt_layer_mua{2})
            for L3=1:length(lkt_layer_mua{3})
                for L4=1:length(lkt_layer_mua{4})
                    lkt_index=find(mua_param_arr(:,1)==lkt_layer_mua{1}(L1) & mua_param_arr(:,2)==lkt_layer_mua{2}(L2) & mua_param_arr(:,3)==lkt_layer_mua{3}(L3) & mua_param_arr(:,4)==lkt_layer_mua{4}(L4));
                    in_place_arr_mua(L1,L2,L3,L4)=lkt_index;
                end
            end
        end
    end
    
    save(fullfile(output_dir,'for_S2.mat'),'in_place_arr','mua_param_arr');
    
    lkt_process_timer=tic;

    %% find the lookup table value for every mus set of the given mua
    temp_mua_param_arr=mua_param_arr; % make the array another shape for better array multiply performance
    temp_mua_param_arr(:,5)=mua_param_arr(:,4)*0.5;
    temp_mua_param_arr=temp_mua_param_arr';
    % divide the mua array into many subarray if it's too large
    num_mua_set=ceil(size(temp_mua_param_arr,2)/max_mua_sameTime); % how many mua subarray
    in_array_mua_param_arr=cell(1,num_mua_set);

    for mua_set_i=1:num_mua_set
        if mua_set_i*max_mua_sameTime<size(temp_mua_param_arr,2)
            this_temp_mua_set=temp_mua_param_arr(:,(mua_set_i-1)*max_mua_sameTime+1:mua_set_i*max_mua_sameTime);
        else
            this_temp_mua_set=temp_mua_param_arr(:,(mua_set_i-1)*max_mua_sameTime+1:end);
        end
        in_array_mua_param_arr{mua_set_i}=reshape(this_temp_mua_set,1,5,size(this_temp_mua_set,2));
    end
    lkt_ref_value_arr_bs=zeros(size(lkt_mus_table,1),size(mua_param_arr,1),num_total); % the the lookup table reflectance value
    lkt_ref_value_arr_as=zeros(size(lkt_mus_table,1),size(mua_param_arr,1),num_total); % the the lookup table reflectance value

    num_SDS=lkt_sim_set.num_SDS;
    
    if exist(fullfile(output_dir,'lkt_ref_value_arr_bs.mat'),'file')
        lkt_ref_value_arr_bs=load(fullfile(output_dir,'lkt_ref_value_arr_bs.mat'));
        lkt_ref_value_arr_bs=lkt_ref_value_arr_bs.lkt_ref_value_arr_bs;
        fprintf('Finish loading lkt_rf_value_arr_bs.mat\n');
    else
        for lkt_index=1:size(lkt_mus_table,1)
            fprintf('Processing lookup mus set %d/%d, Gate: ',lkt_index,size(lkt_mus_table,1));
            temp_PL=load(fullfile(lkt_dir,['sim_' num2str(lkt_index)],'PL_1.mat'));
            for s=1:num_SDS
                SDS_detpt_arr_2(num_gate*(s-1)+1:num_gate*s)=temp_PL.SDS_detpt_arr(:,s)';
            end

            for s=1:num_total
                fprintf(' %d',s);
                temp_lkt_ref_value=zeros(size(mua_param_arr,1),1);
                for mua_set_i=1:num_mua_set
                    temp_detpt_weight=squeeze(exp(-1*sum(SDS_detpt_arr_2{s}(:,1:5).*in_array_mua_param_arr{mua_set_i},2))); % the weight of each photon, not normalized yet
                    if size(SDS_detpt_arr_2{s},1)==1
                        temp_detpt_weight=transpose(temp_detpt_weight/temp_PL.each_photon_weight_arr(ceil(s/10)));
                    else
                        temp_detpt_weight=transpose(sum(temp_detpt_weight,1)/temp_PL.each_photon_weight_arr(ceil(s/10)));
                    end
                    temp_lkt_ref_value(max_mua_sameTime*(mua_set_i-1)+1:max_mua_sameTime*(mua_set_i-1)+length(temp_detpt_weight),1)=temp_detpt_weight; % the reflectance of high NA detector
                end

                lkt_ref_value_arr_bs(lkt_index,:,s)=temp_lkt_ref_value;
            end
            fprintf('\n');
        end
    end
    if ~exist(output_dir,'file')
        mkdir(fullfile(output_dir))
    end
    save(fullfile(output_dir,'lkt_ref_value_arr_bs.mat'),'lkt_ref_value_arr_bs');
    
    lkt_process_timer=toc(lkt_process_timer);
    save(fullfile(output_dir,'lit_process_time.txt'),'lkt_process_timer','-ascii','-tabs');
   
    %% Smooth setting
    close all;
    
    mus_array={'\mu_{s,scalp}','\mu_{s,skull}','\mu_{s,GM}'};
    mua_array={'\mu_{a,scalp}','\mu_{a,skull}','\mu_{a,GM}'};
    tissue_arr={'scalp','skull','GM'};
    SDS_array=[1.5 2.2 2.9 3.6 4.3];
    subplot_arr=[1 4 2 5 3 6];
    
    hs_or_ls=1; % high scattering=0, low scattering=1
    
    linewidth=1.25;
    
    %% Smooth
    
    % for polyfitn
    [x, y, z] = meshgrid(mus_lb(2):25:mus_ub(2),mus_lb(1):25:mus_ub(1), mus_lb(3):25:mus_ub(3));
    x=x(:);
    y=y(:);
    z=z(:);
    
    fig=figure(1); fig.Units='pixels'; fig.Position=[0 0 1920 1080];
    set(fig,'visible','off');
    if exist(fullfile(output_dir,'lkt_ref_value_arr_as.mat'),'file')
        lkt_ref_value_arr_as=load(fullfile(output_dir,'lkt_ref_value_arr_as.mat'));
        lkt_ref_value_arr_as=lkt_ref_value_arr_as.lkt_ref_value_arr_as;
        fprintf('Finish loading lkt_rf_value_arr_as.mat\n');
    else
        for mua_i=1:size(mua_param_arr,1) % choose which mua to fix
            fprintf('Smoothing each mua set %d/%d\n',mua_i,size(mua_param_arr,1));
            for s=1:num_total
                now_SDS=ceil(s/num_gate);
                lkt_value=lkt_ref_value_arr_bs(:,mua_i,s);
                lkt_points_ref=lkt_value(in_place_arr);
    %             lkt_points_ref_arrange(end+1,:,:,:,:)=lkt_points_ref;


                temp_lkt_points_ref=squeeze(lkt_points_ref);

                [m, n, p]=size(temp_lkt_points_ref);
                data_1d=(reshape(temp_lkt_points_ref, [], 1));
                data_1d=-log10(data_1d);
                data_1d(data_1d==inf)=NaN;

                degree=4;
                coefficients=polyfitn([x, y, z], data_1d, degree);
                fitted_values=polyvaln(coefficients, [x, y, z]);
                fitted_values=power(10,-fitted_values);

                fitted_surface=reshape(fitted_values, [m, n, p]);
                smoothed_data(:,:,1,:)=fitted_surface;

                [size1, size2, size3, size4] = size(lkt_points_ref);
                for i=1:length(lkt_value)
                    ind=find(in_place_arr==i);
                    [ind1, ind2, ind3, ind4] = ind2sub([size1, size2, size3, size4], ind);
                    lkt_value_as(1,i)=squeeze(smoothed_data(ind1, ind2, ind3, ind4));
                end

                lkt_ref_value_arr_as(:,mua_i,s)=lkt_value_as;

    %             % Plot as heatmap before and after smooth
    %             figure('Units','pixels','position',[0 0 1920 1080]);
    %             subplot(2,3,1);
    %             cm = [1 0 0;1 1 1; 0 0 1];
    %             cmi = interp1([-100; 0; 100], cm, (-100:100));
    % %             h=heatmap(log10(squeeze(lkt_points_ref(:,:,1,1))),'XData',mus_lb(2):25:mus_ub(2),'YData',mus_lb(1):25:mus_ub(1),'Colormap',cmi ,'CellLabelColor','none','GridVisible','off');
    %             h=heatmap(log10(squeeze(lkt_points_ref(:,:,1,1))),'XData',mus_lb(2):25:mus_ub(2),'YData',mus_lb(1):25:mus_ub(1) ,'CellLabelColor','none','GridVisible','off');
    %             h.NodeChildren(3).YDir='normal';
    %             xlabel('scalp');
    %             ylabel('skull');
    %             title('Reflectance');
    %             
    %             subplot(2,3,2);
    %             h=heatmap(log10(squeeze(lkt_points_ref(:,1,1,:))),'XData',mus_lb(4):25:mus_ub(4),'YData',mus_lb(1):25:mus_ub(1),'CellLabelColor','none','GridVisible','off');
    %             h.NodeChildren(3).YDir='normal';
    %             xlabel('GM');
    %             ylabel('scalp');
    %             title('Reflectance');
    %             
    %             subplot(2,3,3);
    %             h=heatmap(log10(squeeze(lkt_points_ref(1,:,1,:))),'XData',mus_lb(4):25:mus_ub(4),'YData',mus_lb(2):25:mus_ub(2),'CellLabelColor','none','GridVisible','off');
    %             h.NodeChildren(3).YDir='normal';
    %             xlabel('GM');
    %             ylabel('skull');
    %             title('Reflectance');
    %             
    %             subplot(2,3,4);
    %             h=heatmap(log10(squeeze(smoothed_data_1(:,:,1,1))),'XData',mus_lb(2):25:mus_ub(2),'YData',mus_lb(1):25:mus_ub(1) ,'CellLabelColor','none','GridVisible','off');
    %             h.NodeChildren(3).YDir='normal';
    %             xlabel('scalp');
    %             ylabel('skull');
    %             title('Reflectance');
    %             
    %             subplot(2,3,5);
    %             h=heatmap(log10(squeeze(smoothed_data_1(:,1,1,:))),'XData',mus_lb(4):25:mus_ub(4),'YData',mus_lb(1):25:mus_ub(1),'CellLabelColor','none','GridVisible','off');
    %             h.NodeChildren(3).YDir='normal';
    %             xlabel('GM');
    %             ylabel('scalp');
    %             title('Reflectance');
    %             
    %             subplot(2,3,6);
    %             h=heatmap(log10(squeeze(smoothed_data_1(1,:,1,:))),'XData',mus_lb(4):25:mus_ub(4),'YData',mus_lb(2):25:mus_ub(2),'CellLabelColor','none','GridVisible','off');
    %             h.NodeChildren(3).YDir='normal';
    %             xlabel('GM');
    %             ylabel('skull');
    %             title('Reflectance');
    %             
    %             sgtitle(['SDS ' num2str(SDS_array(now_SDS)) ' cm, Gate ' num2str(s-num_gate*(now_SDS-1))]);
    %             if ~exist(fullfile('results','gif'),'file')
    %                 mkdir(fullfile('results','gif'));
    %             end
    %             print(fullfile('results','gif',['compare_SDS' num2str(SDS_array(now_SDS)) 'cm_Gate ' num2str(s-num_gate*(now_SDS-1)) '.png']),'-dpng','-r200');
    %             
    %             % Plot without original
    %             fig=figure('Units','pixels','position',[0 0 1920 540]);
    %             set(fig,'visible','off');
    %             subplot(1,3,1);
    %             h=heatmap(log10(squeeze(smoothed_data_1(:,:,1,1))),'XData',mus_lb(2):25:mus_ub(2),'YData',mus_lb(1):25:mus_ub(1) ,'CellLabelColor','none','GridVisible','off');
    %             h.NodeChildren(3).YDir='normal';
    %             xlabel('scalp');
    %             ylabel('skull');
    %             title('Reflectance');
    %             
    %             subplot(1,3,2);
    %             h=heatmap(log10(squeeze(smoothed_data_1(:,1,1,:))),'XData',mus_lb(4):25:mus_ub(4),'YData',mus_lb(1):25:mus_ub(1),'CellLabelColor','none','GridVisible','off');
    %             h.NodeChildren(3).YDir='normal';
    %             xlabel('GM');
    %             ylabel('scalp');
    %             title('Reflectance');
    %             
    %             subplot(1,3,3);
    %             h=heatmap(log10(squeeze(smoothed_data_1(1,:,1,:))),'XData',mus_lb(4):25:mus_ub(4),'YData',mus_lb(2):25:mus_ub(2),'CellLabelColor','none','GridVisible','off');
    %             h.NodeChildren(3).YDir='normal';
    %             xlabel('GM');
    %             ylabel('skull');
    %             title('Reflectance');
    %             
    %             sgtitle(['SDS ' num2str(SDS_array(now_SDS)) ' cm, Gate ' num2str(s-num_gate*(now_SDS-1))]);
    %             
    %             print(fullfile('results','gif',['SDS' num2str(SDS_array(now_SDS)) 'cm_Gate ' num2str(s-num_gate*(now_SDS-1)) '.png']),'-dpng','-r200');



    %             % Plot as 3d for different SDS and gate
    %             f=figure('Units','pixels','position',[0 0 1920 540]);
    %             set(f,'visible','off');
    %             t=tiledlayout(1,3);
    %             nexttile;
    %             [X,Y] = meshgrid(mus_lb(1):25:mus_ub(1),mus_lb(2):25:mus_ub(2));
    %             surf(X,Y,squeeze(smoothed_data_1(:,:,1,2))','EdgeColor','None');
    %             xlabel('scalp');
    %             ylabel('skull');
    %             zlabel('Reflectance');
    %             
    %             nexttile;
    %             [X,Y] = meshgrid(mus_lb(1):25:mus_ub(1),mus_lb(4):25:mus_ub(4));
    %             surf(X,Y,squeeze(smoothed_data_1(:,2,1,:))','EdgeColor','None');
    %             xlabel('scalp');
    %             ylabel('gray matter');
    %             zlabel('Reflectance');
    %             
    %             nexttile;
    %             [X,Y] = meshgrid(mus_lb(2):25:mus_ub(2),mus_lb(4):25:mus_ub(4));
    %             surf(X,Y,squeeze(smoothed_data_1(2,:,1,:))','EdgeColor','None');
    %             xlabel('skull');
    %             ylabel('gray matter');
    %             zlabel('Reflectance');
    %             
    %             title(t,['SDS ' num2str(SDS_array(now_SDS)) ' cm, Gate ' num2str(s-num_gate*(now_SDS-1))]);
    %             print(fullfile('results',['SDS' num2str(SDS_array(now_SDS)) 'cm_Gate ' num2str(s-num_gate*(now_SDS-1)) '.png']),'-dpng','-r200');


    %             % Prepare to plot
    %             lkt_points_ref_smooth(end+1,:,:,:,:)=smoothed_data;
    % %             lkt_points_ref_smooth_1(end+1,:,:,:,:)=smoothed_data_1;
    % %             lkt_points_ref_smooth_2(end+1,:,:,:,:)=smoothed_data_2;
    %             
    %             if rem(s,num_gate)==0
    %                 now_SDS=s/num_gate;
    % 
    %                 for mus=[1 2 4]
    %                     % smallest mus
    %                     if mus==1
    %                         if hs_or_ls==0
    %                             to_plot_orig=squeeze(lkt_points_ref_arrange(:,:,8,1,10))';
    %                             to_plot_smooth=squeeze(lkt_points_ref_smooth(:,:,8,1,10))';
    %                         elseif hs_or_ls==1
    %                             to_plot_orig=squeeze(lkt_points_ref_arrange(:,:,2,1,2))';
    %                             to_plot_smooth=squeeze(lkt_points_ref_smooth(:,:,2,1,2))';
    %                         end
    %                     elseif mus==2
    %                         if hs_or_ls==0
    %                             to_plot_orig=squeeze(lkt_points_ref_arrange(:,7,:,1,10))';
    %                             to_plot_smooth=squeeze(lkt_points_ref_smooth(:,7,:,1,10))';
    % %                             to_plot_true_sk=dtof_arrange_sk(:,9+(now_SDS-1)*num_gate:8+now_SDS*num_gate);
    %                         elseif hs_or_ls==1
    %                             to_plot_orig=squeeze(lkt_points_ref_arrange(:,2,:,1,2))';
    %                             to_plot_smooth=squeeze(lkt_points_ref_smooth(:,2,:,1,2))';
    % %                             to_plot_true_sk=dtof_arrange_sk(:,9+(now_SDS-1)*num_gate:8+now_SDS*num_gate);
    %                         end
    %                     elseif mus==4
    %                         if hs_or_ls==0
    %                             to_plot_orig=squeeze(lkt_points_ref_arrange(:,7,8,1,:))';
    %                             to_plot_smooth=squeeze(lkt_points_ref_smooth(:,7,8,1,:))';
    % %                             to_plot_true_gm=dtof_arrange_gm(:,9+(now_SDS-1)*num_gate:8+now_SDS*num_gate);
    %                         elseif hs_or_ls==1
    %                             to_plot_orig=squeeze(lkt_points_ref_arrange(:,2,2,1,:))';
    %                             to_plot_smooth=squeeze(lkt_points_ref_smooth(:,2,2,1,:))';
    % %                             to_plot_true_gm=dtof_arrange_gm(:,9+(now_SDS-1)*num_gate:8+now_SDS*num_gate);
    %                         end
    %                     end
    %                     
    %                     % Plot comparison of original and smooth data in each SDS and gate
    %                     fig2=figure('Units','pixels','position',[0 0 1920 1080]);
    % %                     set(fig2,'visible','off');
    %                     ti=tiledlayout(2,5);
    %                     for g=1:num_gate
    %                         nexttile;
    %                         plot(mus_lb(mus):25:mus_ub(mus),to_plot_orig(:,g),'--o','Linewidth',linewidth);
    %                         hold on
    %                         plot(mus_lb(mus):25:mus_ub(mus),to_plot_smooth(:,g),'-o','Linewidth',linewidth);
    %                         hold on
    % %                         plot(mus_lb(mus):25:mus_ub(mus),to_plot_smooth_1(:,g),'-o');
    % %                         hold on
    % %                         if mus==4
    % %                             plot(mus_lb(mus):25:mus_ub(mus),to_plot_true_gm(:,g),'-o','Linewidth',linewidth);
    % %                             hold on
    % %                         elseif mus==2
    % %                             plot(mus_lb(mus):25:mus_ub(mus),to_plot_true_sk(:,g),'-o','Linewidth',linewidth);
    % %                             hold on
    % %                         end
    %                         x1=mus_lb(mus)+25;
    %                         x2=mus_ub(mus)-25;
    %                         yLimits = ylim;
    %                         xPatch = [x1, x2, x2, x1];
    %                         yPatch = [min(ylim), min(ylim),max(ylim), max(ylim)];
    %                         p=patch(xPatch, yPatch, [0.69, 0.93, 0.93],'FaceAlpha',0.3,'EdgeColor','none');
    %                         p.Annotation.LegendInformation.IconDisplayStyle = 'off';
    %                         uistack(p,"bottom");    
    %                         xlabel(mus_array{mus});
    %                         ylabel('reflectance');
    %                         ylim(yLimits);
    %                         title(['Gate' num2str(g)]);
    %                         legend('1E10','polynomial','2E11','Location','northwest');
    %                     end
    %                     title(ti,['SDS ' num2str(SDS_array(now_SDS)) ' cm']);
    %                     
    %                     if ~exist(fullfile('smooth_result',subject_name),'file')
    %                         mkdir(fullfile('smooth_result',subject_name))
    %                     end
    %                     
    %                     if hs_or_ls==0
    %                         print(fig2,fullfile('smooth_result',subject_name,['hs_mus_' num2str(mus) '_SDS ' num2str(SDS_array(now_SDS)) 'cm_smooth_result.png']),'-dpng','-r200');
    %                     elseif hs_or_ls==1
    %                         print(fig2,fullfile('smooth_result',subject_name,['ls_mus_' num2str(mus) '_SDS ' num2str(SDS_array(now_SDS)) 'cm_smooth_result.png']),'-dpng','-r200');
    %                     end
    %                 end
    %                 lkt_points_ref_arrange=[];
    %                 lkt_points_ref_smooth=[];
    %                 lkt_points_ref_smooth_1=[];
    %                 lkt_points_ref_smooth_2=[];
    %             end
            end
        end
        save(fullfile(output_dir,'lkt_ref_value_arr_as.mat'),'lkt_ref_value_arr_as');
    end
    
    %% Plot for changing mus, fixed mua (1D)
    hs_or_ls_or_mid=0; % high scattering=0, low scattering=1, medium scattering=2
    
    mua_i=size(lkt_ref_value_arr_bs,2);
    % load true answer
    dtof_arrange_gm=[];
    mus_table = load(fullfile(subject_name,'mus_table.txt'));
    for sim = 1:size(mus_table,1)
        dtof=load(fullfile(subject_name,['DTOF_' num2str(sim) '.mat']));
        dtof_arrange_gm(end+1,:)=dtof.to_save;
    end
    
    dtof_arrange_ls=dtof_arrange_gm(1:11,:);
    dtof_arrange_hs=dtof_arrange_gm(12:22,:);
    dtof_arrange_ms=dtof_arrange_gm(1:11,:);
    
    % plot the simulation result separately

    for s=21:50
        if rem(s,num_gate)==1
            figure('Units','pixels','position',[0 0 1920 1080]);
            ti=tiledlayout(3,10);
        end

        now_SDS=ceil(s/num_gate);
        now_gate=s-(now_SDS-1)*num_gate;

        lkt_value=lkt_ref_value_arr_bs(:,mua_i,s);
        lkt_4D_bs=lkt_value(in_place_arr);

        lkt_value=lkt_ref_value_arr_as(:,mua_i,s);
        lkt_4D_as=lkt_value(in_place_arr);

        if rem(s,num_gate)<11 % only plot gate1~6 rem(s,num_gate)>0 && rem(s,num_gate)<7
            for mus=1:3
                if mus==1
                    if hs_or_ls_or_mid==0
                        to_plot_orig=squeeze(lkt_4D_bs(:,8,1,10));
                        to_plot_smooth=squeeze(lkt_4D_as(:,8,1,10));
                    elseif hs_or_ls_or_mid==1
                        to_plot_orig=squeeze(lkt_4D_bs(:,2,1,2))';
                        to_plot_smooth=squeeze(lkt_4D_as(:,2,1,2))';
                    elseif hs_or_ls_or_mid==2
                        to_plot_orig=squeeze(lkt_4D_bs(:,5,1,6))';
                        to_plot_smooth=squeeze(lkt_4D_as(:,5,1,6))';
                    end
                elseif mus==2
                    if hs_or_ls_or_mid==0
                        to_plot_orig=squeeze(lkt_4D_bs(7,:,1,10))';
                        to_plot_smooth=squeeze(lkt_4D_as(7,:,1,10))';
    %                             to_plot_true_sk=dtof_arrange_sk(:,9+(now_SDS-1)*num_gate:8+now_SDS*num_gate);
                    elseif hs_or_ls_or_mid==1
                        to_plot_orig=squeeze(lkt_4D_bs(2,:,1,2))';
                        to_plot_smooth=squeeze(lkt_4D_as(2,:,1,2))';
    %                             to_plot_true_sk=dtof_arrange_sk(:,9+(now_SDS-1)*num_gate:8+now_SDS*num_gate);
                     elseif hs_or_ls_or_mid==2
                        to_plot_orig=squeeze(lkt_4D_bs(4,:,1,6))';
                        to_plot_smooth=squeeze(lkt_4D_as(4,:,1,6))';
    %                             to_plot_true_sk=dtof_arrange_sk(:,9+(now_SDS-1)*num_gate:8+now_SDS*num_gate);
                    end
                elseif mus==3
                    if hs_or_ls_or_mid==0
                        to_plot_orig=squeeze(lkt_4D_bs(7,8,1,:))';
                        to_plot_smooth=squeeze(lkt_4D_as(7,8,1,:))';
                        to_plot_true_gm=dtof_arrange_hs(:,8+now_gate+(now_SDS-1)*num_gate);
                    elseif hs_or_ls_or_mid==1
                        to_plot_orig=squeeze(lkt_4D_bs(2,2,1,:))';
                        to_plot_smooth=squeeze(lkt_4D_as(2,2,1,:))';
                        to_plot_true_gm=dtof_arrange_ls(:,8+now_gate+(now_SDS-1)*num_gate);
                    elseif hs_or_ls_or_mid==2
                        to_plot_orig=squeeze(lkt_4D_bs(4,5,1,:))';
                        to_plot_smooth=squeeze(lkt_4D_as(4,5,1,:))';
                        to_plot_true_gm=dtof_arrange_ms(:,8+now_gate+(now_SDS-1)*num_gate);
                    end
                end
                nexttile(now_gate+10*(mus-1));
                plot(mus_lb(mus):25:mus_ub(mus),to_plot_orig,'--o','Linewidth',linewidth);
                hold on
                plot(mus_lb(mus):25:mus_ub(mus),to_plot_smooth,'-o','Linewidth',linewidth);
                if mus==3
                    hold on
                    plot(mus_lb(mus):25:mus_ub(mus),to_plot_true_gm,'-o','Linewidth',linewidth);
                end

                % plot underground with real boundary
                x1=mus_lb(mus)+25;
                x2=mus_ub(mus)-25;
                yLimits = ylim;
                xPatch = [x1, x2, x2, x1];
                yPatch = [min(ylim), min(ylim),max(ylim), max(ylim)];
                p=patch(xPatch, yPatch, [0.69, 0.93, 0.93],'FaceAlpha',0.3,'EdgeColor','none');
                p.Annotation.LegendInformation.IconDisplayStyle = 'off';
                uistack(p,"bottom");

                xlabel(mus_array{mus});
                ylabel('reflectance');
                ylim(yLimits);
                title(['Gate' num2str(now_gate)]);
%                 legend('1E10','polynomial','2E11','Location','northwest');
                title(ti,['SDS ' num2str(SDS_array(now_SDS)) ' cm']);
            end
        end


        if rem(s,num_gate)==0
            lgd=legend('1E10','polynomial','2E11','Orientation','horizontal');
            lgd.Layout.Tile='south';
            if hs_or_ls_or_mid==0
               print(fullfile(output_dir,['hs_SDS_' num2str(SDS_array(now_SDS)) 'cm_smooth_result.png']),'-dpng','-r200');
            elseif hs_or_ls_or_mid==1
               print(fullfile(output_dir,['ls_SDS_' num2str(SDS_array(now_SDS)) 'cm_smooth_result.png']),'-dpng','-r200');
            elseif hs_or_ls_or_mid==2
               print(fullfile(output_dir,['ms_SDS_' num2str(SDS_array(now_SDS)) 'cm_smooth_result.png']),'-dpng','-r200');
            end
        end
    end
    
    
    
    %% Plot for changing mus, fixed mua (2D)
    % Plot as 3d for specific SDS and gate (please choose one 's' to plot, or you will get many figures for each 's')
    plot_fig=0;
    
    if plot_fig==1
        for s=35
            now_SDS=ceil(s/num_gate);
            now_gate=s-(now_SDS-1)*num_gate;

            lkt_value=lkt_ref_value_arr_bs(:,mua_i,s);
            lkt_4D_bs=lkt_value(in_place_arr);

            lkt_value=lkt_ref_value_arr_as(:,mua_i,s);
            lkt_4D_as=lkt_value(in_place_arr);

            % before smoothing
            [X,Y] = meshgrid(mus_lb(1):25:mus_ub(1),mus_lb(2):25:mus_ub(2));
            zmax=max(lkt_4D_bs(:));
            zmin=min(lkt_4D_bs(:));
            for i=1:size(lkt_4D_bs,4)
                f=figure('Units','pixels','position',[0 0 640 540]);
                surf(X,Y,squeeze(lkt_4D_bs(:,:,1,i))','EdgeColor','None');
                title(['GM = ' num2str(mus_lb(3)+25*(i-1)) ' cm^{-1}']);
                xlabel('scalp');
                ylabel('skull');
                zlabel('Reflectance');
                zlim([zmin zmax]);
                print(fullfile(output_dir,['b_fixed_GM_' num2str(i) '_SDS' num2str(SDS_array(now_SDS)) 'cm_Gate ' num2str(s-num_gate*(now_SDS-1)) '.png']),'-dpng','-r200');
            end

            [X,Y] = meshgrid(mus_lb(1):25:mus_ub(1),mus_lb(3):25:mus_ub(3));
            for i=1:size(lkt_4D_bs,2)
                f=figure('Units','pixels','position',[0 0 640 540]);
                surf(X,Y,squeeze(lkt_4D_bs(:,i,1,:))','EdgeColor','None');
                title(['skull = ' num2str(mus_lb(2)+25*(i-1)) ' cm^{-1}']);
                xlabel('scalp');
                ylabel('gray matter');
                zlabel('Reflectance');
                zlim([zmin zmax]);
                print(fullfile(output_dir,['b_fixed_skull_' num2str(i) '_SDS' num2str(SDS_array(now_SDS)) 'cm_Gate ' num2str(s-num_gate*(now_SDS-1)) '.png']),'-dpng','-r200');
            end

            [X,Y] = meshgrid(mus_lb(2):25:mus_ub(2),mus_lb(3):25:mus_ub(3));
            for i=1:size(lkt_4D_bs,1)
                f=figure('Units','pixels','position',[0 0 640 540]);
                surf(X,Y,squeeze(lkt_4D_bs(i,:,1,:))','EdgeColor','None');
                title(['scalp = ' num2str(mus_lb(1)+25*(i-1)) ' cm^{-1}']);
                xlabel('skull');
                ylabel('gray matter');
                zlabel('Reflectance');
                zlim([zmin zmax]);
                print(fullfile(output_dir,['b_fixed_scalp_' num2str(i) '_SDS' num2str(SDS_array(now_SDS)) 'cm_Gate ' num2str(s-num_gate*(now_SDS-1)) '.png']),'-dpng','-r200');
            end

            % after smoothing
            [X,Y] = meshgrid(mus_lb(1):25:mus_ub(1),mus_lb(2):25:mus_ub(2));
            zmax=max(lkt_4D_as(:));
            zmin=min(lkt_4D_as(:));
            for i=1:size(lkt_4D_as,4)
                f=figure('Units','pixels','position',[0 0 640 540]);
                surf(X,Y,squeeze(lkt_4D_as(:,:,1,i))','EdgeColor','None');
                title(['GM = ' num2str(mus_lb(3)+25*(i-1)) ' cm^{-1}']);
                xlabel('scalp');
                ylabel('skull');
                zlabel('Reflectance');
                zlim([zmin zmax]);
                print(fullfile(output_dir,['a_fixed_GM_' num2str(i) '_SDS' num2str(SDS_array(now_SDS)) 'cm_Gate ' num2str(s-num_gate*(now_SDS-1)) '.png']),'-dpng','-r200');
            end

            [X,Y] = meshgrid(mus_lb(1):25:mus_ub(1),mus_lb(3):25:mus_ub(3));
            for i=1:size(lkt_4D_as,2)
                f=figure('Units','pixels','position',[0 0 640 540]);
                surf(X,Y,squeeze(lkt_4D_as(:,i,1,:))','EdgeColor','None');
                title(['skull = ' num2str(mus_lb(2)+25*(i-1)) ' cm^{-1}']);
                xlabel('scalp');
                ylabel('gray matter');
                zlabel('Reflectance');
                zlim([zmin zmax]);
                print(fullfile(output_dir,['a_fixed_skull_' num2str(i) '_SDS' num2str(SDS_array(now_SDS)) 'cm_Gate ' num2str(s-num_gate*(now_SDS-1)) '.png']),'-dpng','-r200');
            end

            [X,Y] = meshgrid(mus_lb(2):25:mus_ub(2),mus_lb(3):25:mus_ub(3));
            for i=1:size(lkt_4D_as,1)
                f=figure('Units','pixels','position',[0 0 640 540]);
                surf(X,Y,squeeze(lkt_4D_as(i,:,1,:))','EdgeColor','None');
                title(['scalp = ' num2str(mus_lb(1)+25*(i-1)) ' cm^{-1}']);
                xlabel('skull');
                ylabel('gray matter');
                zlabel('Reflectance');
                zlim([zmin zmax]);
                print(fullfile(output_dir,['a_fixed_scalp_' num2str(i) '_SDS' num2str(SDS_array(now_SDS)) 'cm_Gate ' num2str(s-num_gate*(now_SDS-1)) '.png']),'-dpng','-r200');
            end
        end
    end

    
    
    %% truth for changing mua, fixed mus
%     true_dir=fullfile(lookup_table_arr,'KB_high_scattering');
%     for lkt_index=11 %size(lkt_mus_table,1)
%         fprintf('Processing lookup mus set %d/%d, Gate: ',lkt_index,size(lkt_mus_table,1));
%         temp_PL=load(fullfile(true_dir,['sim_' num2str(lkt_index)],'PL_1.mat'));
%         for s=1:num_SDS
%             SDS_detpt_arr_2(num_gate*(s-1)+1:num_gate*s)=temp_PL.SDS_detpt_arr(:,s)';
%         end
% 
%         for s=1:num_total
%             fprintf(' %d',s);
%             temp_lkt_ref_value=zeros(size(mua_param_arr,1),1);
%             for mua_set_i=1:num_mua_set
%                 temp_detpt_weight=squeeze(exp(-1*sum(SDS_detpt_arr_2{s}(:,1:5).*in_array_mua_param_arr{mua_set_i},2))); % the weight of each photon, not normalized yet
%                 temp_detpt_weight=transpose(sum(temp_detpt_weight,1)/temp_PL.each_photon_weight_arr(ceil(s/10)));
%                 temp_lkt_ref_value(max_mua_sameTime*(mua_set_i-1)+1:max_mua_sameTime*(mua_set_i-1)+length(temp_detpt_weight),1)=temp_detpt_weight; % the reflectance of high NA detector
%             end
%             temp_lkt_ref_value=temp_lkt_ref_value(in_place_arr_mua);
%             lkt_ref_value_true_changing_mua(s,:,:,:,:)=temp_lkt_ref_value;
%         end
%         
%         fprintf('\n');
%     end
%     
%     %% For changing mua, fixed mus
%     close all;
%     
%     lkt_points_ref_arrange=[];
%     lkt_points_ref_smooth=[];
%     lkt_points_ref_smooth_1=[];
%     lkt_points_ref_smooth_2=[];
%     
%     orig=[];
%     smooth=[];
%     smooth_1=[];
%     
%     % for polyfitn
%     [x, y, z] = meshgrid(mus_lb(2):25:mus_ub(2),mus_lb(1):25:mus_ub(1), mus_lb(4):25:mus_ub(4));
%     x=x(:);
%     y=y(:);
%     z=z(:);
%     
%     for s=21:50 %21:50 %1:num_total
%         for mua_i=1:size(mua_param_arr,1)
%             lkt_value=lkt_ref_value_arr(:,mua_i,s);
%             lkt_points_ref=lkt_value(in_place_arr);
%             orig(end+1)=lkt_points_ref(7,8,1,10); % choose which mus to fix
%             
%             % smooth
%             temp_lkt_points_ref=squeeze(lkt_points_ref);
%            
%             [m, n, p]=size(temp_lkt_points_ref);
%             data_1d=(reshape(temp_lkt_points_ref, [], 1));
%             data_1d=-log10(data_1d);
%             data_1d(data_1d==inf)=NaN;
%             
%             degree=4;
%             coefficients=polyfitn([x, y, z], data_1d, degree);
%             fitted_values=polyvaln(coefficients, [x, y, z]);
%             fitted_values=power(10,-fitted_values);
%             
%             fitted_surface=reshape(fitted_values, [m, n, p]);
%             smoothed_data(:,:,1,:)=fitted_surface;
%             
% %             degree=6;
% %             coefficients=polyfitn([x, y, z], data_1d, degree);
% %             fitted_values=polyvaln(coefficients, [x, y, z]);
% %             fitted_values=power(10,-fitted_values);
% %             
% %             fitted_surface=reshape(fitted_values, [m, n, p]);
% %             smoothed_data_1(:,:,1,:)=fitted_surface;
%             
%             
%             smooth(end+1)=smoothed_data(7,8,1,10);
% %             smooth_1(end+1)=smoothed_data_1(7,8,1,10);
% 
%         end
%         
%         lkt_points_ref_arrange(end+1,:,:,:,:)=orig(in_place_arr_mua);
%         lkt_points_ref_smooth(end+1,:,:,:,:)=smooth(in_place_arr_mua);
% %         lkt_points_ref_smooth_1(end+1,:,:,:,:)=smooth_1(in_place_arr_mua);
%         
%         orig=[];
%         smooth=[];
%         smooth_1=[];
%         
%         if rem(s,10)==0
%             now_SDS=s/num_gate;
%             for mua=[1 2 4]
%                 figure('Units','pixels','position',[0 0 1920 1080]);
%                 ti=tiledlayout(2,5);
%                 % smallest mus
%                 if mua==1
%                     to_plot_orig=squeeze(lkt_points_ref_arrange(:,:,1,1,1))';
%                     to_plot_smooth=squeeze(lkt_points_ref_smooth(:,:,1,1,1))';
%                     to_plot_true=squeeze(lkt_ref_value_true_changing_mua(s-9:s,:,1,1,1))';
%                 elseif mua==2
%                     to_plot_orig=squeeze(lkt_points_ref_arrange(:,1,:,1,1))';
%                     to_plot_smooth=squeeze(lkt_points_ref_smooth(:,1,:,1,1))';
%                     to_plot_true=squeeze(lkt_ref_value_true_changing_mua(s-9:s,1,:,1,1))';
%                 elseif mua==4
%                     to_plot_orig=squeeze(lkt_points_ref_arrange(:,1,1,1,:))';
%                     to_plot_smooth=squeeze(lkt_points_ref_smooth(:,1,1,1,:))';
%     %                 to_plot_smooth_1=squeeze(lkt_points_ref_smooth_1(:,1,1,1,:))';
%     %                 to_plot_smooth_2=squeeze(lkt_points_ref_smooth_2(:,1,1,1,:))';
%                     to_plot_true=squeeze(lkt_ref_value_true_changing_mua(s-9:s,1,1,1,:))';
%                 end
% 
%                 % Plot comparison of original and smooth data in each SDS and gate
%                 for g=1:num_gate
%                     nexttile;
%                     plot(mua_lb(mua):0.05:mua_ub(mua),to_plot_orig(:,g),'--o','LineWidth',linewidth);
%                     hold on
%                     plot(mua_lb(mua):0.05:mua_ub(mua),to_plot_smooth(:,g),'-o','LineWidth',linewidth);
%                     hold on 
%                     plot(mua_lb(mua):0.05:mua_ub(mua),to_plot_true(:,g),'-o','LineWidth',linewidth);
%                     
% %                     xLimits=xlim;
% %                     yLimits=ylim;
% %                     xPatch=[mua_lb(mua), mua_ub(mua), mua_ub(mua), mua_lb(mua)];
% %                     yPatch=[min(ylim), min(ylim),max(ylim), max(ylim)];
% %                     p=patch(xPatch, yPatch, [0.69, 0.93, 0.93],'FaceAlpha',0.3,'EdgeColor','none');
% %                     p.Annotation.LegendInformation.IconDisplayStyle = 'off';
% %                     uistack(p,"bottom");
%                     
%                     xlabel(mua_array{mua});
%                     ylabel('reflectance');
%                     title(['Gate' num2str(g)]);
% %                     xlim(xLimits);
% %                     ylim(yLimits);
%                     legend('1E10','degree=4','2E11','Location','northwest');
%                 end
%                 title(ti,['SDS ' num2str(SDS_array(now_SDS)) ' cm']);
%                 print(fullfile('results',['hs_mua_' num2str(mua) '_SDS ' num2str(SDS_array(now_SDS)) 'cm_smooth_result.png']),'-dpng','-r200');
%             end
%             lkt_points_ref_arrange=[];
%             lkt_points_ref_smooth=[];
%             lkt_points_ref_smooth_1=[];
%             lkt_points_ref_smooth_2=[];
%         end
%     end
end

disp('Done!');

