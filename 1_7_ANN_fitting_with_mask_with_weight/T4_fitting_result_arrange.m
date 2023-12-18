%{
Arrange the fitting result of the random generated target spectrum

Ting-Yi Kuo
Last update: 2023/08/21
%}

clc;clear;close all; clearvars -global;

global lambda cw_net cw_param_range tr_net tr_param_range 

%% param
input_dir='test_fitting_2023-12-12-11-59-43'; % please move the fitting folders into this folder first.
subject_name_arr={'KB','WH','ZJ'}; %,'WH','WW'
num_anser_to_generate=10; % number of target spec (true answer)
num_error_to_generate=1; % number of adding noise to the same, the first one will have no error
times_to_fitting=20; % number of fitting using different init value

fitting_dir='fitting_tr1234_gate12'; % the fitting folder

do_plot_anyPlot=0; % =1 to plot any plot
do_plot_individual_fitting=1; % =1 to plot the result of each individual fitting

% CW setting
num_SDS_cw=6; % how many SDS are in the target spectrum
SDS_sim_correspond_cw=[1 2 3 4 5 6]; % the SDS index in the simulated spectrum corresponding to each SDS in the target spectrum, SDS_sim_correspond(x)=y means 'measure x = sim y'

% TR setting
num_SDS_tr=5;
num_gate=10;
SDS_sim_correspond_tr=[1 2 3 4 5]; % the SDS index in the simulated spectrum corresponding to each SDS in the target spectrum, SDS_sim_correspond(x)=y means 'measure x = sim y'


% other setting
num_fitted_param=13; % the number of the fitted parameters

SDS_length_arr_cw=[0.8 1.5 2.12 3 3.35 4.5]; % cm
SDS_length_arr_tr=[1.5 2.2 2.9 3.6 4.3];

fitting_wl_file='fitting_wl.txt'; % the file in the epsilon folder containing the fitting wavelength
fitting_wl_tr=810;

model_dir='model_arrange'; % the folder containing the arranged ANN file

fontSize=12;
lineWidth=2;
lgdFontSize=12;
lgdNumCol=7;

subplot_height=250; % pixel, the height of subplot
subplot_width=450; % pixel, the width of subplot
left_spacing=80; % pixel, the space between subplot and the left things
right_spacing=50; % pixel, the space right of the last column of subplot
upper_spacing=70; % pixel, the space between subplot and the upper things
lower_spacing=30; % pixel, the space below the legend box
legend_height=50; % pixel, the height of legend box

plot_n_col=1; % the column number of subplot
plot_n_row=1; % the row number of subplot

%% init
fitting_wl=load(fullfile('epsilon',fitting_wl_file));
fitting_wl=[fitting_wl;fitting_wl_tr];
lambda=unique(fitting_wl);


mkdir(input_dir,'arrangement');

% load the param and OP answer
param_answer_arr=load(fullfile(input_dir,'answers','param_answer_arr.txt'));
OP_answer_arr=[];
for i=1:num_anser_to_generate
    temp_OP_arr=load(fullfile(input_dir,'answers',['OP_ans_' num2str(i) '.txt']));
    OP_answer_arr(:,:,i)=interp1(temp_OP_arr(:,1),temp_OP_arr(:,2:end),lambda);
end

%% main

for sbj=1:length(subject_name_arr)
    % make output folder and find fitting folder
    mkdir(fullfile(input_dir,'arrangement',subject_name_arr{sbj}));
    
    % load ANN model
    load(fullfile(model_dir,[subject_name_arr{sbj} '_cw_model.mat'])); % cw_net, cw_param_range
    load(fullfile(model_dir,[subject_name_arr{sbj} '_tr_model.mat'])); % tr_net, tr_param_range

    for target_i=1:num_anser_to_generate
        for error_i=1:num_error_to_generate
            fprintf('Process %s target %d, error %d, fitting ',subject_name_arr{sbj},target_i,error_i);

            temp_target_folder=fullfile(input_dir,subject_name_arr{sbj},['target_' num2str(target_i) '_' num2str(error_i)],fitting_dir);
            
            % skip the processed folder
            if exist(fullfile(temp_target_folder,'fitting_results.csv'),'file')
                fprintf('\n');
                continue;
            end

            target_spec=load(fullfile(input_dir,subject_name_arr{sbj},['target_' num2str(target_i) '_' num2str(error_i)],'target_spec.txt'));
            target_dtof=load(fullfile(input_dir,subject_name_arr{sbj},['target_' num2str(target_i) '_' num2str(error_i)],'target_dtof.txt'));

            init_param_arr=load(fullfile(temp_target_folder,'init_arr.txt'));
            load(fullfile(temp_target_folder,'fitting_info.mat'));
            if cw_flag && ~tr_flag
                SDS_choosed_tr=[];
            elseif ~cw_flag && tr_flag
                SDS_choosed_cw=[];
            end
            
            fitting_result=[];
            fitting_op_arr=[];
            fitted_ann_spec_arr_=[];

            fun_init_param_to_mu_spec(); % load the epsilon for fitting wl

            interped_target_spec=interp1(target_spec(:,1),target_spec(:,2:end),lambda); % for calculate error
            
            fitted_ann_spec_arr=[];
            init_ann_spec_arr=[];
            calculated_param_error_arr=[];
            
            for j=1:times_to_fitting
                fprintf('%d ',j);
                
                % loaf the fitting result
                temp_route=load(fullfile(temp_target_folder,['fitting_' num2str(j)],'fitRoute.txt'));
                fitting_result(j,1:num_fitted_param+num_SDS_cw+num_SDS_tr+3)=temp_route(end,:); % the fitted param
                fitting_op_arr(:,:,j)=load(fullfile(temp_target_folder,['fitting_' num2str(j)],'fitted_mu.txt'));
                
                %% forward the init param and calculate the error
                [OP_arr,~]=fun_param_to_mu(init_param_arr(j,1:num_fitted_param),0);
                
                % CW
                if cw_flag
                    init_ann_spec=fun_ANN_forward(OP_arr,0);
                    init_ann_spec_arr(:,:,j)=init_ann_spec(:,SDS_sim_correspond_cw);

                    calculated_init_error=sqrt(mean((init_ann_spec_arr(:,:,j)./interped_target_spec-1).^2,1));
                    init_error(1)=sqrt(mean(calculated_init_error(SDS_choosed_cw).^2));
                else
                    calculated_init_error=NaN*ones(1,num_SDS_cw);
                    init_error(1)=NaN;
                end
                
                % TR
                if tr_flag
                    temp_OP_arr=[lambda OP_arr];
                    OP_arr=interp1(temp_OP_arr(:,1),temp_OP_arr(:,2:end),fitting_wl_tr);
                    init_ann_dtof_1D=fun_ANN_forward(OP_arr,1);
                    for i=1:num_SDS_tr
                        init_ann_dtof(:,i)=init_ann_dtof_1D((i-1)*num_gate+1:i*num_gate)';
                    end
                    init_ann_dtof_arr(:,:,j)=init_ann_dtof(:,SDS_sim_correspond_tr);

                    select_gate=sum(mask,1);
                    calculated_init_error_tr=init_ann_dtof_arr(:,:,j)./target_dtof-1;
                    calculated_init_error_tr=calculated_init_error_tr.*mask;  % choose gate
                    calculated_init_error_tr=sqrt(sum(calculated_init_error_tr.^2,1)./select_gate);
                    init_error(2)=sqrt(mean(calculated_init_error_tr(SDS_choosed_tr).^2));
                else
                    calculated_init_error_tr=NaN*ones(1,num_SDS_tr);
                    init_error(2)=NaN;
                end
                
                % all
                temp_init_error=init_error;
                temp_init_error(isnan(temp_init_error))=0;
                init_error_all=(temp_init_error(1).*weight(1)+temp_init_error(2).*weight(2))/sum(weight);

                init_param_arr(j,num_fitted_param+1:num_fitted_param+num_SDS_cw+num_SDS_tr+3)=[calculated_init_error calculated_init_error_tr init_error init_error_all];
                
                
                %% Calculate the fitted error
                if cw_flag
                    fitted_ann_spec_arr(:,:,j)=load(fullfile(temp_target_folder,['fitting_' num2str(j)],'fitted_spec.txt'));
                    calculated_fitted_error=sqrt(mean((fitted_ann_spec_arr(:,:,j)./interped_target_spec-1).^2,1));
                    fitted_error(1)=sqrt(mean(calculated_fitted_error(SDS_choosed_cw).^2));
                else
                    % use ANN to get supposed spectra
                    fitted_op=squeeze(fitting_op_arr(:,:,j));
                    fitted_ann_spec=fun_ANN_forward(fitted_op,0);
                    fitted_ann_spec_arr(:,:,j)=fitted_ann_spec(:,SDS_sim_correspond_cw);
                    calculated_fitted_error=sqrt(mean((fitted_ann_spec_arr(:,:,j)./interped_target_spec-1).^2,1));
                    fitted_error(1)=NaN;
                end

                if tr_flag
                    fitted_ann_dtof_arr_2D(:,:,j)=load(fullfile(temp_target_folder,['fitting_' num2str(j)],'fitted_dtof.txt'));
                    calculated_fitted_error_tr=fitted_ann_dtof_arr_2D(:,:,j)./target_dtof(:,:)-1;
                    calculated_fitted_error_2D_arr(:,:,j)=calculated_fitted_error_tr; % save errors
                    calculated_fitted_error_tr=calculated_fitted_error_tr.*mask;  % choose gate
                    calculated_fitted_error_tr=sqrt((sum(calculated_fitted_error_tr.^2,1)./select_gate));
                    fitted_error(2)=sqrt(mean(calculated_fitted_error_tr(SDS_choosed_tr).^2));
                else
                    % use ANN to get supposed DTOFs
                    fitted_op=squeeze(fitting_op_arr(:,:,j));
                    fitted_op=interp1(lambda,fitted_op,fitting_wl_tr);
                    fitted_ann_dtof=fun_ANN_forward(fitted_op,1);
                    for i=1:num_SDS_tr
                        fitted_ann_dtof_arr(:,i)=fitted_ann_dtof((i-1)*num_gate+1:(i-1)*num_gate+num_gate)';
                    end
                    fitted_ann_dtof_arr_2D(:,:,j)=fitted_ann_dtof_arr;
                    
                    % Calculate for drawing 
                    calculated_fitted_error_tr=fitted_ann_dtof_arr_2D(:,:,j)./target_dtof(:,:)-1;
                    calculated_fitted_error_2D_arr(:,:,j)=calculated_fitted_error_tr; % save errors
                    calculated_fitted_error_tr=sqrt((mean(calculated_fitted_error_tr.^2,1)));
                    fitted_error(2)=NaN;
                end
                
                temp_fitted_error=fitted_error;
                temp_fitted_error(isnan(temp_fitted_error))=0;
                fitted_error_all=(temp_fitted_error(1).*weight(1)+temp_fitted_error(2).*weight(2))/sum(weight);
                    

                %% save the errors
                fitting_result(j,num_fitted_param+1:num_fitted_param+num_SDS_cw+num_SDS_tr+3)=[calculated_fitted_error calculated_fitted_error_tr fitted_error fitted_error_all];
                init_param_arr(j,num_fitted_param+1:num_fitted_param+num_SDS_cw+num_SDS_tr+3)=[calculated_init_error calculated_init_error_tr init_error init_error_all];
                
                % calculate the param and OP answer
                if ~cw_flag && tr_flag % fitting only TR
                    fitting_op=interp1(lambda,fitting_op_arr(:,:,j),fitting_wl_tr,'pchip');
                    OP_answer=interp1(lambda,OP_answer_arr(:,:,target_i),fitting_wl_tr,'pchip');
                    calculated_OP_error=sqrt(mean((fitting_op./OP_answer-1).^2,1));
                else
                    calculated_OP_error=sqrt(mean((fitting_op_arr(:,:,j)./OP_answer_arr(:,:,target_i)-1).^2,1));
                end
                    
                calculated_OP_error=calculated_OP_error([1 2 3 4 7 8]); % only calculate L1 2 4
                calculated_OP_error_save=[calculated_OP_error sqrt(mean(calculated_OP_error.^2))];
                fitting_result(j,num_fitted_param+num_SDS_cw+num_SDS_tr+4:num_fitted_param+num_SDS_cw+num_SDS_tr+3+size(calculated_OP_error_save,2))=calculated_OP_error_save;
            end
            fprintf('\n');


            %% output CSV
            % sort according to reflectance error
            [~,error_index]=sort(fitting_result(:,num_fitted_param+num_SDS_cw+num_SDS_tr+3),'ascend');
            fitting_result_sort=[error_index fitting_result(error_index,:)];
            init_param_arr_sort=init_param_arr(error_index,:);
%             calculated_param_error_arr_sort=calculated_param_error_arr(error_index,:);
            save(fullfile(temp_target_folder,'fitting_result_sort.txt'),'fitting_result_sort','-ascii','-tabs');
            save(fullfile(temp_target_folder,'init_param_arr_sort.txt'),'init_param_arr_sort','-ascii','-tabs');
%             save(fullfile(temp_target_folder,'calculated_param_error_arr_sort.txt'),'calculated_param_error_arr_sort','-ascii','-tabs');
            
            fid=fopen(fullfile(temp_target_folder,'fitting_results.csv'),'w');
            
            % header
            % for fitted param
            fprintf(fid,',A1,K1,A2,K2,A4,K4,hc1,sto2_1,hc2,sto2_2,hc4,sto2_4,mel');
            % for SDS error
            for s=1:num_SDS_cw
                fprintf(fid,',error %d',s);
            end
            
            for s=1:num_SDS_tr
                fprintf(fid,',error-%d',s);
            end
            
            fprintf(fid,',cw_choosed_error,tr_choosed_error,error_choosed');
            % for fitted OP error
            for L=[1 2 4]
                fprintf(fid,',mua %d,mus %d',L,L);
            end
            fprintf(fid,'\n');
            
            % print the true answer
            fprintf(fid,'answer');
            fprintf(fid,',%f',param_answer_arr(target_i,:));
            fprintf(fid,'\n');

            % fitting result
            for j=error_index'
                fprintf(fid,'Fitting %d',j);
                % for fitted param
                fprintf(fid,',%f',fitting_result(j,1:num_fitted_param));
                % for SDS error
                fprintf(fid,',%.2f%%',fitting_result(j,num_fitted_param+1:num_fitted_param+num_SDS_cw+num_SDS_tr+3)*100);
                % for OP error
                fprintf(fid,',%.2f%%',fitting_result(j,num_fitted_param+num_SDS_cw+num_SDS_tr+4:end)*100);
                fprintf(fid,'\n');
                % for param error
%                 fprintf(fid,'param error');
%                 fprintf(fid,',%.2f%%',calculated_param_error_arr(j,:)*100);
%                 fprintf(fid,'\n');
            end
            fclose(fid);



            if do_plot_anyPlot
                %% plot the OPs
                fprintf('Plot fitted OP\n');
                legend_arr={'answer'};
                for j=1:times_to_fitting
                    legend_arr{end+1}=['rank ' num2str(times_to_fitting+1-j) ', ' num2str(fitting_result(error_index(times_to_fitting+1-j),end)*100,'%.2f%%')];
                end

                fig=figure('Units','pixels','position',[0 0 1920 1080]);
                ti=tiledlayout('flow','TileSpacing','Compact','Padding','Compact');
    %             set(fig,'visible','off');
                op_name={'\mu_{a','\mu_{s'};
                layer_name_arr={'scalp','skull','CSF','GM'};
                colormap_arr=jet(times_to_fitting);
                colormap_arr=colormap_arr(end:-1:1,:);
                for L=[1 2 4]
                    for opi=1:2
                        nexttile;
                        hold on;
                        plot(lambda,OP_answer_arr(:,2*(L-1)+opi,target_i),'LineWidth',lineWidth);
                        for j=times_to_fitting:-1:1
                            plot(lambda,fitting_op_arr(:,2*(L-1)+opi,error_index(j)),'--','Color',colormap_arr(j,:),'LineWidth',1.5); %,'LineWidth',lineWidth
                        end
                        % scale to ub and lb
                        ylim([cw_param_range(2,(opi-1)*4+L) cw_param_range(1,(opi-1)*4+L)]);
                        xlabel('wavelength(nm)');
                        ylabel([op_name{opi} ',' layer_name_arr{L} '}(1/cm)']);
                        set(gca,'fontsize',fontSize, 'FontName', 'Times New Roman');
                        grid on;
                    end
                end
                
                lgd=legend(legend_arr,'fontsize',lgdFontSize);
                lgd.NumColumns=lgdNumCol;
                lgd.Layout.Tile = 'south';

                title(ti,subject_name_arr{sbj});

                print(fullfile(temp_target_folder,'fitted_mu.png'),'-dpng','-r200');
    %                 close all;
                
                %% plot the fitted spec together (for CW)
                fprintf('Plot the fitted spec\n');
                legend_arr={'target'};
                for j=16:20 % 1:times_to_fitting
                    legend_arr{end+1}=['rank ' num2str(times_to_fitting+1-j) ', ' num2str(fitting_result(error_index(times_to_fitting+1-j),num_fitted_param+num_SDS_cw+num_SDS_tr+1)*100,'%.2f%%')];
                end

                fig=figure('Units','pixels','position',[0 0 1980 1080]);
                ti=tiledlayout('flow','TileSpacing','Compact','Padding','Compact');
%                 set(fig,'visible','off');
                op_name={'mua','mus'};
                colormap_arr=jet(5); % times_to_fitting
                colormap_arr=colormap_arr(end:-1:1,:);

                for s=1:num_SDS_cw
                    nexttile;
                    plot(target_spec(:,1),target_spec(:,s+1),'LineWidth',lineWidth);
                    hold on;
                    for j=5:-1:1 %times_to_fitting
                       plot(lambda,fitted_ann_spec_arr(:,s,error_index(j)),'--','Color',colormap_arr(j,:),'LineWidth',lineWidth);
                        hold on;
                    end
                    title(['SDS ' num2str(SDS_length_arr_cw(s)) ' cm']);
                    ylabel('reflectance');
                    set(gca,'fontsize',fontSize, 'FontName', 'Times New Roman');
    %                 set(gca,'Unit','pixels','Position',[0 0 480 400]);
                    grid on;
                end

                lgd=legend(legend_arr,'fontsize',lgdFontSize);
                lgd.NumColumns=lgdNumCol;
                lgd.Layout.Tile = 'south';
                
                title(ti,subject_name_arr{sbj});
   
                print(fullfile(temp_target_folder,'fitted_spec.png'),'-dpng','-r200');
%                 close all;
                
                %% plot the fitted dtof error together (for TR)
                fprintf('Plot the fitted error\n');
                legend_arr={'target'};
                for j=16:20 %1:times_to_fitting
                    legend_arr{end+1}=['rank ' num2str(times_to_fitting+1-j) ', ' num2str(fitting_result(error_index(times_to_fitting+1-j),num_fitted_param+num_SDS_cw+num_SDS_tr+2)*100,'%.2f%%')];
                end
                
                fig=figure('Units','pixels','position',[0 0 1920 1080]);
                ti=tiledlayout('flow','TileSpacing','Compact','Padding','Compact');
%                 set(fig,'visible','off');
                op_name={'mua','mus'};
                colormap_arr=jet(5);%times_to_fitting
                colormap_arr=colormap_arr(end:-1:1,:);
                

                for s=1:num_SDS_tr
                    nexttile;
                    r=[1:0.1:num_gate];
                    y = zeros(length(r),1);
                    plot(r,y,'LineWidth',lineWidth);
                    hold on;
                    for j=5:-1:1 %times_to_fitting
                        plot(1:1:num_gate,100*calculated_fitted_error_2D_arr(:,s,error_index(j)),'--','Color',colormap_arr(j,:),'LineWidth',lineWidth);
                        hold on;
                    end
    %                 ylim([-10 10]);
                    title(['SDS ' num2str(SDS_length_arr_tr(s)) ' cm']);
                    ylabel('error(%)');
                    xlim([1 10]);
                    xticks(1:1:10);
                    set(gca,'fontsize',fontSize, 'FontName', 'Times New Roman');
    %                 set(gca,'Unit','pixels','Position',[0 0 480 400]);
                    grid on;
                    
                    % plot underground with fitting range
                    if ismember(s,SDS_choosed_tr)
                        first=0;
                        last=0;
                        for g=1:num_gate
                            if mask(g,s)==1
                                first=g;
                                break
                            end
                        end

                        for g=num_gate:-1:1
                            if mask(g,s)==1
                                last=g;
                                break
                            end
                        end

                        x1=first;
                        x2=last;
                        yLimits = ylim;
                        xPatch = [x1, x2, x2, x1];
                        yPatch = [min(ylim), min(ylim),max(ylim), max(ylim)];
                        p=patch(xPatch, yPatch, [0.69, 0.93, 0.93],'FaceAlpha',0.4,'EdgeColor','none');
                        p.Annotation.LegendInformation.IconDisplayStyle = 'off';
                        uistack(p,"bottom");
                        ylim(yLimits);
                    end
                end

                lgd=legend(legend_arr,'fontsize',lgdFontSize);
                lgd.NumColumns=lgdNumCol;
                lgd.Layout.Tile = 'south';
                
                title(ti,subject_name_arr{sbj});

                print(fullfile(temp_target_folder,'fitted_DTOF_error.png'),'-dpng','-r200');
%                 close all;

                % plot the real error of each SDS and gate with heatmap
                figure;
                h=heatmap(round(100*sqrt(squeeze(calculated_fitted_error_2D_arr(:,:,error_index(1)).^2)),2));
%                 h.YDisplayData=flipud(h.YDisplayData);
                xlabel('SDS');
                ylabel('gate');
                title('RMSPE')
                annotation('textarrow',[1,1],[0.5,0.5],'string','(%)', ...
                'HeadStyle','none','LineStyle','none','HorizontalAlignment','center','TextRotation',90);
                print(fullfile(temp_target_folder,'fitted_DTOF_error_detail.png'),'-dpng','-r200');
                
            end
            
            %% plot the fitted error and the init error, also the OP error and the fitted error
            if do_plot_anyPlot
                fprintf('Plot the fitted param\n');
                fig=figure('Units','pixels','position',[0 0 1600 1200]);
                set(fig,'visible','off');
                ti=tiledlayout('flow','TileSpacing','compact','Padding','none');
                % fitted error and init error
                nexttile();
                plot(init_param_arr(:,end),fitting_result(:,num_fitted_param+num_SDS_cw+num_SDS_tr+3),'o');
                grid on;
                xlabel('init error');
                ylabel('fitted error');
                xxlim=xlim(); xxlim(1)=0; xlim(xxlim);
                yylim=ylim(); yylim(1)=0; ylim(yylim);
                set(gca,'fontsize',fontSize, 'FontName', 'Times New Roman');
                % error for each param
%                 param_name_arr={'hc_1','sto2_1','hc_2','sto2_2','hc_4','sto2_4','mel'};
%                 for param_i=7:num_fitted_param % only plot the mua param
%                     nexttile();
%                     plot(fitting_result(:,num_fitted_param+num_SDS+1),calculated_param_error_arr(:,param_i),'o');
%                     grid on;
%                     xlabel('spectral error');
%                     ylabel('percent error');
%                     yline(0,'b','LineWidth',lineWidth);
%                     title(param_name_arr{param_i-6});
%                     set(gca,'fontsize',fontSize, 'FontName', 'Times New Roman');
%                 end
                % OP error
                param_name_arr={'\mu_{a,scalp}','\mu_{s,scalp}','\mu_{a,skull}','\mu_{s,skull}','\mu_{a,GM}','\mu_{s,GM}','all OPs'};
                for param_i=1:7
                    nexttile();
                    plot(fitting_result(:,num_fitted_param+num_SDS_cw+num_SDS_tr+3),fitting_result(:,num_fitted_param+num_SDS_cw+num_SDS_tr+3+param_i),'o');
                    grid on;
                    xlabel('spectral error');
                    ylabel('RMSPE');
                    title(param_name_arr{param_i});
                    yylim=ylim(); yylim(1)=0; ylim(yylim);
                    set(gca,'fontsize',fontSize, 'FontName', 'Times New Roman');
                end
                print(fullfile(temp_target_folder,'fitted_param.png'),'-dpng','-r200');
%                 close all;
            end

            %% plot the fitted spec
            if do_plot_individual_fitting && do_plot_anyPlot
                
                if cw_flag
    %                 for jj=1:times_to_fitting % the rank
                    for jj=1:1 % the rank
                        fprintf('Plot %s fitting %d\n',subject_name_arr{sbj},jj);
                        j=error_index(jj); % the index of fitting
                        fig=figure('Units','pixels','position',[0 0 1600 900]);
                        set(fig,'visible','off');
                        ti=tiledlayout('flow','TileSpacing','compact','Padding','none');

                        for s=1:num_SDS_cw
                            nexttile();
                            hold on;
                            plot(target_spec(:,1),target_spec(:,s+1),'LineWidth',lineWidth);
                            plot(lambda,init_ann_spec_arr(:,s,j),'LineWidth',lineWidth);
                            plot(lambda,fitted_ann_spec_arr(:,s,j),'LineWidth',lineWidth);
                            grid on;

                            lgd=legend({'target',['init, err=' num2str(init_param_arr(j,num_fitted_param+s)*100,'%.2f%%')],['fitted, err=' num2str(fitting_result(j,num_fitted_param+s)*100,'%.2f%%')]},'Location','southoutside');
                            lgd.NumColumns=3;
                            title(['SDS ' num2str(SDS_length_arr_cw(s)) ' cm']);
                            set(gca,'fontsize',fontSize, 'FontName', 'Times New Roman');
                        end

                        title(ti,[strrep(subject_name_arr{sbj},'_',' ') ' target ' num2str(target_i) ' error ' num2str(error_i) ', rank ' num2str(jj) ' = fitting ' num2str(j) ', fitting error = ' num2str(fitting_result(j,num_fitted_param+num_SDS_cw+num_SDS_tr+1)*100,'%.2f%%')], 'FontName', 'Times New Roman');
                        print(fullfile(temp_target_folder,['fitted_spec_r' num2str(jj) '_f' num2str(j) '.png']),'-dpng','-r200');
    %                     close all;
                    end
                end
                
                if tr_flag
                    %         for jj=1:times_to_fitting % the rank
                    for jj=1:1 % the rank
                        fprintf('Plot %s fitting %d\n',subject_name_arr{sbj},jj);
                        j=error_index(jj); % the index of fitting
                        fig=figure('Units','pixels','position',[0 0 1600 900]);
                        set(fig,'visible','off');
                        ti=tiledlayout('flow','TileSpacing','compact','Padding','none');

                        for s=1:num_SDS_tr
                            nexttile();
                            semilogy(1:1:num_gate,target_dtof(:,s),'LineWidth',lineWidth);
                            hold on;
                            semilogy(1:1:num_gate,init_ann_dtof_arr(:,s,j),'LineWidth',lineWidth);
                            hold on;
                            semilogy(1:1:num_gate,fitted_ann_dtof_arr_2D(:,s,j),'LineWidth',lineWidth);
                            grid on;

                            lgd=legend({'target',['init, err=' num2str(init_param_arr(j,num_fitted_param+num_SDS_cw+s)*100,'%.2f%%')],['fitted, err=' num2str(fitting_result(j,num_fitted_param+num_SDS_cw+s)*100,'%.2f%%')]},'Location','southoutside');
                            lgd.NumColumns=3;
        %                         title(['SDS ' num2str(SDS_length_arr(s)) ' cm']);
                            set(gca,'fontsize',fontSize, 'FontName', 'Times New Roman');
                        end

                        title(ti,[strrep(subject_name_arr{sbj},'_',' ') ' target ' num2str(target_i) ' error ' num2str(error_i) ', rank ' num2str(jj) ' = fitting ' num2str(j) ', fitting error = ' num2str(fitting_result(j,num_fitted_param+num_SDS_cw+num_SDS_tr+2)*100,'%.2f%%')], 'FontName', 'Times New Roman');
                        print(fullfile(temp_target_folder,['fitted_spec_r' num2str(jj) '_f' num2str(j) '_tr.png']),'-dpng','-r200');
    %                     close all;
                    end
                end
            end
        end
    end
end

disp('Done!');