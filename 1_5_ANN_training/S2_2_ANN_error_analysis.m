%{
Find the ANN error and the relationship with the inputs

Benjamin Kao
Last update: 2020/12/23
%}

clc;clear;close all;

%% param
subject_name={'KB'};
input_dir_arr={'KB_2023-10-11-20-35-47'};
max_sequence_length=250;
num_SDS=5;
num_gate=10;

SDS_dist_arr=[1.5 2.2 2.9 3.6 4.3]; % cm
num_error_to_check=400;

fontSize=14;
marker_size=20;
line_width=2;

for sbj_i=1:length(input_dir_arr)
    input_dir=input_dir_arr{sbj_i};
    %% init
    load(fullfile(input_dir,'ANN__train_info.mat'));
    param_range=load(fullfile(input_dir,'param_range.txt'));
    param_range=param_range(:,[1 2 4 5 6 8]);

    true_input_arr=zeros(max_sequence_length*length(testing_input),size(testing_input{1},1));
    for i=1:length(testing_input)
        true_input_arr((i-1)*max_sequence_length+1:i*max_sequence_length,:)=testing_input{i}';
    end

    true_input_arr=normalize_param(true_input_arr,param_range,2);
    
    load(fullfile('../1_3_MCX_lookup_table/results/',['no_photons_SDS45_gate1_' subject_name{sbj_i} '.mat'])); % SDS_to_examine, gate_to_examine, record_table
    
    %% main
    for s=1:num_SDS
        for g=1:num_gate
            ind=g+(s-1)*num_gate;
            fprintf('Analysis %s SDS %d gate %d\n',input_dir,s,g);
            fig=figure('Units','pixels','position',[0 0 2000 600]);
            set(fig, 'visible', 'off');
    %         ti=tiledlayout('flow','TileSpacing','compact','Padding','none');
            ti=tiledlayout(1,4,'TileSpacing','compact','Padding','none');

            [SDS_error,error_index]=sort(abs(error(:,ind)),'descend');
            to_check_index=error_index(1:num_error_to_check);
            to_check_error=error(to_check_index,ind);
            error_input=true_input_arr(to_check_index,:);
            error_answer=testing_output(to_check_index,ind); % the answer in the training data
            error_output=reflectance_arr(to_check_index,ind); % the output of the ANN

            % find if the params are close to the boundary
            % examine both mus and mua
            input_closeTo_ub=abs((error_input./param_range(1,:)-1))<0.02;
            input_closeTo_lb=abs((error_input./param_range(2,:)-1))<0.02;
            input_closeTo_boundary=sum(input_closeTo_ub+input_closeTo_lb,2);

            nexttile();
            plot(abs(to_check_error),input_closeTo_boundary,'.','MarkerSize',marker_size);
            xlabel('error');
            ylabel('number of parameter close to boundary');
            ylim([0 length(param_range)]);
            set(gca,'fontsize',fontSize, 'FontName', 'Times New Roman');

            % find the relationship with the original answer
            nexttile();
            plot(error_answer,error_output,'.','MarkerSize',marker_size);
            hold on;
            plot([0 max(error_answer)],[0 max(error_answer)],'LineWidth',line_width);
            xlabel('true answer');
            ylabel('ANN output');
            set(gca,'fontsize',fontSize, 'FontName', 'Times New Roman');

            % find the percentage of the true answer in the training data
            [~,all_reflectance_index]=sort(testing_output(:,ind),'ascend');
            in_reflectance_percentage=[];
            for i=1:num_error_to_check
                in_reflectance_percentage(i,1)=find(all_reflectance_index==to_check_index(i));
            end
            in_reflectance_percentage=in_reflectance_percentage/size(reflectance_arr,1);
            nexttile();
            plot(in_reflectance_percentage,'.','MarkerSize',marker_size);
            xlabel('error rank');
            ylabel('true answer percentage in all training data');
            set(gca,'fontsize',fontSize, 'FontName', 'Times New Roman');

            title(ti,['SDS = ' num2str(SDS_dist_arr(s)) ' cm Gate = ' num2str(g) ' , largest ' num2str(num_error_to_check) ' error'],'fontsize',fontSize, 'FontName', 'Times New Roman');
            drawnow;
            print(fullfile(input_dir,['error_analysis_SDS_' num2str(SDS_dist_arr(s)) 'cm_gate_' num2str(g) '.png']),'-dpng','-r200');
            
            % Analyze no photon effect on error
%             if ismember(s,SDS_to_examine) && ismember(g,gate_to_examine)
                exists_arr=false(size(error_input, 1), 1);
                
                if isempty(record_table{g,s})
                    true_count=0;
                    false_count=400;
                else
                    for i = 1:size(error_input, 1)
                        exists_arr(i) = any(ismember(error_input(i,4:6),record_table{g,s}(:,[1 2 4]),'rows'));
                    end
                    true_count = sum(exists_arr);
                    false_count = numel(exists_arr)-true_count;
                end

                counts = [true_count, false_count];
                labels = {'True', 'False'};
                nexttile;
                ax=gca();
                p=pie(counts, labels);
                newColors = [...
                    0.3010, 0.7450, 0.9330
                    0, 0.4470, 0.7410;];
                ax.Colormap = newColors;
                title('Receive no photons');

                print(fullfile(input_dir,['error_analysis_SDS_' num2str(SDS_dist_arr(s)) 'cm_gate_' num2str(g) '_pie.png']),'-dpng','-r200');
%             end
        end
 
        close all;
    end
    
    
    
end

disp('Done!');

% normalize the parameters to [0,1]
% direction: if =1, normalize the input; if =2, denormalize the input
function output=normalize_param(input,param_range,direction)
    param_scaling=param_range(1,:)-param_range(2,:);
    if direction==1
        output=(input-param_range(2,:))./param_scaling;
    elseif direction==2
        output=input.*param_scaling+param_range(2,:);
    else
        assert(false,'Function ''normalize_spec'' param ''direction'' Error!');
    end
end