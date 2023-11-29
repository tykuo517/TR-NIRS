%{
Train the ANN of head model using the lookup table predicted value.

Benjamin Kao 
Last update: 2020/12/22
Edit by Ting-Yi Kuo 2023/10/12
%}

clc;clear;close all;

%% param

subject_name_arr={'WH'}; % the name of the subject
input_mother_dir='../1_4_MCX_lookup_table_dataGen'; % the folder of the lookup table forward training data
input_dir_arr={'WH_2023-11-14-21-53-01'}; % the lut folder of each subject

for sbj_i=1:length(subject_name_arr)
    clearvars -except subject_name_arr input_mother_dir input_dir_arr sbj_i;
    close all;
    
    subject_name=subject_name_arr{sbj_i};
    input_dir=fullfile(input_mother_dir,input_dir_arr{sbj_i});
    num_SDS=5;
    num_gate=10;
    do_normalize=1; % if =1, do normalization to spec and param
    max_sequence_length=250; % how many data sets in one batch, to prevent use too much GPU memory 250

    testing_rate=0.15; % the ratio of testing data
    validation_rate=0.1; % the ratio of validation data

    %% load
    fprintf('loading training data\n');
    all_data=load(fullfile(input_dir,'all_param_arr.mat'));
    all_data=all_data.all_param_arr;
    param_range=load(fullfile(input_dir,'param_range.txt'));

    %% init
    output_dir=[subject_name '_' datestr(datetime('now'),'yyyy-mm-dd-HH-MM-ss')];
    mkdir(output_dir);

    rand_seed=now;
    setdemorandstream(rand_seed);
    save(fullfile(output_dir,['rand_seed.txt']),'rand_seed','-ascii');
    
    % shuffle
    all_data=all_data(randperm(size(all_data,1)),:);

    % delete nan and minus value
    positive_index=find(~sum(((all_data(:,9:end)<=0) + (isnan(all_data(:,9:end)))) ,2)>0);
    all_data=all_data(positive_index,:);

    % all_data=all_data(1:floor(size(all_data,1)/20),:);
    all_data=all_data(1:floor(size(all_data,1)/max_sequence_length)*max_sequence_length,:);


    if do_normalize==1
        input_param=normalize_param(all_data(:,1:8),param_range,1);
        save(fullfile(output_dir,'param_range.txt'),'param_range','-ascii','-tabs');
    %     input_param=normalize_param_2(all_data(:,1:8),param_range,1);
        [output_spec,spec_range]=normalize_spec(all_data(:,9:end),1);
        save(fullfile(output_dir,'spec_range.mat'),'spec_range');
    else
        input_param=all_data(:,1:8);
        output_spec=all_data(:,9:end);
    end

    % backup this script
    copyfile([mfilename '.m'],fullfile(output_dir,[mfilename '.m']));
    
    
    input_param = input_param(:,[1 2 4 5 6 8]);

    %% training

    % divide the data into many mini batch, especially the validation data and testing data
    num_training_group=floor((size(input_param,1)*(1-testing_rate-validation_rate))/max_sequence_length);
    num_validation_group=floor(num_training_group/(1-testing_rate-validation_rate)*validation_rate);
    num_testing_group=floor(num_training_group/(1-testing_rate-validation_rate)*testing_rate);
    assert(max_sequence_length*(num_training_group+num_validation_group+num_testing_group)<=size(input_param,1),'ERROR while doing the divide');

%     training_data_input=input_param(1:num_testing_group*max_sequence_length,:);
%     training_data_output=output_spec(1:num_testing_group*max_sequence_length,:);
    training_data_input=input_param(1:num_training_group*max_sequence_length,:);
    training_data_output=output_spec(1:num_training_group*max_sequence_length,:);
    validation_data=cell(1,2);
    validation_data{1,1}=input_param(num_training_group*max_sequence_length+1:(num_validation_group+num_training_group)*max_sequence_length,:)';
    validation_data{1,2}=output_spec(num_training_group*max_sequence_length+1:(num_validation_group+num_training_group)*max_sequence_length,:)';
    testing_input=cell(num_testing_group,1);
    for i=1:num_testing_group
        testing_input{i,1}=input_param((i+num_training_group+num_validation_group-1)*max_sequence_length+1:(i+num_training_group+num_validation_group)*max_sequence_length,:)';
    end
    testing_output=all_data((num_training_group+num_validation_group)*max_sequence_length+1:(num_training_group+num_validation_group+num_testing_group)*max_sequence_length,9:end);


    %     'InitialLearnRate',1e-8, ...

    options = trainingOptions('adam', ...
        'Shuffle','every-epoch', ...
        'ExecutionEnvironment','gpu', ...
        'InitialLearnRate',1e-5, ...
        'LearnRateSchedule','piecewise', ...
        'LearnRateDropFactor',0.95, ...
        'LearnRateDropPeriod',10, ...
        'ValidationData',validation_data, ...
        'SequenceLength',max_sequence_length, ...
        'MaxEpochs',300, ...
        'ValidationFrequency',num_training_group, ... %num_training_group
        'ValidationPatience',50, ... %20
        'MiniBatchSize',128, ... %128
        'VerboseFrequency',1, ...
        'Plots','training-progress');
%         'SequenceLength',max_sequence_length, ...

    
    layers=[sequenceInputLayer(6,'Name','')
%         fullyConnectedLayer(650)
        fullyConnectedLayer(850) %850
        leakyReluLayer
    %     reluLayer
%         fullyConnectedLayer(450)
        fullyConnectedLayer(550) %550
        leakyReluLayer
    %     reluLayer
        fullyConnectedLayer(300) %300
        leakyReluLayer
    %     reluLayer
        fullyConnectedLayer(150) %150
        leakyReluLayer
    %     reluLayer
        fullyConnectedLayer(50)
        regressionLayer];

%     gpuDevice(1);
    [net,info]=trainNetwork(training_data_input',training_data_output',layers,options);

    save(fullfile(output_dir,['ANN_model.mat']),'net');

    reflectance_arr=zeros(num_testing_group*max_sequence_length,num_SDS*num_gate);
    for i=1:num_testing_group
        reflectance=predict(net,testing_input{i,1});
        reflectance_arr((i-1)*max_sequence_length+1:i*max_sequence_length,:)=reflectance';
    end
    if do_normalize==1
        reflectance_arr=normalize_spec(reflectance_arr,2,spec_range);
    end
    error=double(reflectance_arr)./double(testing_output)-1;
    error_mean=mean(abs(error),1);
    error_std=std(error,[],1);
    error_rmspe=sqrt(mean(error.^2,1));
    
    save(fullfile(output_dir,'ANN__train_info.mat'),'info','error','testing_input','testing_output','reflectance_arr');

    for s=1:num_SDS
        figure('Units','pixels','position',[0 0 1920 1080]);
        ti=tiledlayout('flow','TileSpacing','compact','Padding','none');
        for g=1:num_gate
            nexttile();
            hist(error(:,10*(s-1)+g),50);
            xlabel('Error');
            ylabel('testing data number');
            title({['Gate ' num2str(g) ' testing error histogram'],['mean error = ' num2str(error_mean(10*(s-1)+g)*100,'%.2f%%') ', std error = ' num2str(error_std(10*(s-1)+g)*100,'%.2f%%') ', rmspe = ' num2str(error_rmspe(10*(s-1)+g)*100,'%.2f%%')]});

            fprintf('---- error= %2.4f%%, rmse= %2.4f%%, std=%2.4f%%\n',error_mean(10*(s-1)+g)*100,error_rmspe(10*(s-1)+g)*100,error_std(10*(s-1)+g)*100);
        end
        fprintf('\n--------------------------\n\n')
        title(ti,[subject_name ' SDS' num2str(s) ' training result']);
        saveas(gcf,fullfile(output_dir,[subject_name '_SDS' num2str(s) '_testing_error_histogram.png']));
    end
    
    error_record=[error_mean' error_std' error_rmspe'];

    save(fullfile(output_dir,'error_record.txt'),'error_record','-ascii','-tabs');
    % close all;
end

disp('Done!');

%% functions

% direction: if =1, normalize the input; if =2, denormalize the input
function [output,spec_range]=normalize_spec(input,direction,spec_range)
    if direction==1
        output=-log10(input);
        spec_range=[];
%         max_value=max(output);
%         min_value=min(output);
%         spec_range=[max_value;min_value];
%         output=(output-min_value)./(max_value-min_value);
    elseif direction==2
%         input=input.*(spec_range(1,:)-spec_range(2,:))+spec_range(2,:);
        output=power(10,-input);
        spec_range=[];
    else
        assert(false,'Function ''normalize_spec'' param ''direction'' Error!');
    end
end

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

% normalize the parameters to [0,1]
% direction: if =1, normalize the input; if =2, denormalize the input
function output=normalize_param_2(input,param_range,direction)
    param_scaling=param_range(1,:)-param_range(2,:);
    param_mean=mean(param_range,1);
    if direction==1
        output=(input-param_mean)./param_scaling*2;
    elseif direction==2
        output=input.*param_scaling/2+param_mean;
    else
        assert(false,'Function ''normalize_spec'' param ''direction'' Error!');
    end
end