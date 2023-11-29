%{
Use the smooth mcx lookup table to generate training data for ANN.

This code only generate mus, which is different from original mus set, to
get more data through interpolation.
For generate mua please refer to 'S1_smooth_lookup_table.m'.

Please run 'S1_smooth_lookup_table.m' before running this code.


Ting-Yi Kuo
Last update: 2023/10/16
%}

clc;clear;close all;

%% param
lookup_table_arr='/home/md703/Documents/Ty/TR-NIRS/1_3_MCX_lookup_table'; % the dir containing the unmerged lookup table
subject_name_arr={'WH'}; % the name of the subjects

num_layer=4; % number of layer to random
setting_r=0.2; % the radius of true detector, in mm

% mus boundary around 800nm  
mus_ub=[225 200 23 250]; % 1/cm
mus_lb=[100 50 23 50]; % 1/cm

% actual simulation of mus boundary
% mus_ub=[250 225 23 275]; % 1/cm, skip CSF
% mus_lb=[75 25 23 25]; % 1/cm, skip CSF

mua_ub=[0.45 0.3 0.042 0.4]; % 1/cm
mua_lb=[0.1 0.1 0.042 0.1]; % 1/cm


% about random
num_random=1000; % how many number to random, for mua.  For mus, it's this number + number of original lookup table set
normal_cutoff=2; % if the random value is not in [-a a], than re-generate a random number

num_SDS=5;
num_gate=10;
num_total=num_SDS*num_gate;

max_mua_sameTime=10; % how many mua set to calculate at the same time, use smaller value for smaller memory consumption

test_mode=1; % =0 to generate the whole training data; =1 or more to generate the result of lookup table using the testing parameters


for sbj_i=1:length(subject_name_arr)
    
    subject_name=subject_name_arr{sbj_i};
    input_dir=fullfile('results_smooth',subject_name);
    
    load(fullfile(input_dir,'lkt_ref_value_arr_as.mat'));
    load(fullfile(input_dir,'for_S2.mat')); % in_place_arr, mua_param_arr
    mua_param_arr=mua_param_arr(:,1:4);

    %% init
    lkt_dir=fullfile(lookup_table_arr,subject_name); % the dir containing the unmerged lookup table
    
    % make the output dir name
    if test_mode==0
        output_dir=[subject_name '_' datestr(datetime('now'),'yyyy-mm-dd-HH-MM-ss')];
    else
        output_dir=[subject_name '_test' num2str(test_mode) '_' datestr(datetime('now'),'yyyy-mm-dd-HH-MM-ss')];
    end
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
    save(fullfile(output_dir,'param_range.txt'),'to_save','-ascii','-tabs');

    if test_mode==0
        %% random generate the mus array, and combine with the original lookup table mus
        random_index_arr=normrnd(0,1,num_random,num_layer);
        while sum(random_index_arr(:)>normal_cutoff)~=0 && sum(random_index_arr(:)<-normal_cutoff)~=0
            random_index_arr(random_index_arr>normal_cutoff)=normrnd(0,1,length(find(random_index_arr>normal_cutoff)),1);
            random_index_arr(random_index_arr<-normal_cutoff)=normrnd(0,1,length(find(random_index_arr<-normal_cutoff)),1);
        end

        mus_param_arr=zeros(size(random_index_arr));
        for L=1:num_layer
            mus_param_arr(:,L)=(random_index_arr(:,L)+normal_cutoff)/(2*normal_cutoff)*(mus_ub(L)-mus_lb(L))+mus_lb(L);
        end
        % add the known correct mus points in actual boundary
        index=[];
        for j=1:num_layer
            index(:,j)=lkt_mus_table(:,j)<mus_lb(:,j) | lkt_mus_table(:,j)>mus_ub(:,j); % either parameters out of actual boundary 
        end
        index=sum(index,2);
        temp_lkt_mus_table=lkt_mus_table(~index,:);
        mus_param_arr(end+1:end+size(temp_lkt_mus_table,1),:)=temp_lkt_mus_table(:,1:num_layer); % add the known correct mus points
    end

    %% test
    if test_mode==1
        layer_mus={100,50,23,[50:10:250]};
        mus_param_arr=[];
        for i=1:length(layer_mus{1})
            for j=1:length(layer_mus{4})
                mus_param_arr=[mus_param_arr; layer_mus{1}(i) layer_mus{2}(i) layer_mus{3}(i) layer_mus{4}(j)];
            end
        end
        
        % specify the index of mua combination you want to use in mua_param_arr
        index=1331;
        mua_param_arr=mua_param_arr(index,:); %[0.45 0.3 0.042 0.4]
        mua_param_arr=repmat(mua_param_arr,size(mus_param_arr,1),1);
        
    elseif test_mode==2
        op=load('OP_sim_sen.txt');
        mus_param_arr=op(:,5:8);
        mua_param_arr=op(:,1:4);
    elseif test_mode==3
        op=load('OPs_to_sim_11/toSim_OP_66.txt');
        mus_param_arr=op(:,2:2:8);
        mua_param_arr=op(:,1:2:8);
    end


    %% interpolate the lookup table
    mua_param_arr=mua_param_arr(:,1:4);
    if test_mode==0
        all_param_arr=zeros(size(mua_param_arr,1)*size(mus_param_arr,1),2*num_layer+num_total);
    elseif test_mode>=1
        all_param_arr=zeros(size(mua_param_arr,1),2*num_layer+num_total);
    end

    for mua_i=1:size(mua_param_arr,1)
        fprintf('Calculating mua %d/%d, Gate ',mua_i,size(mua_param_arr,1));
        for s=1:num_total
            fprintf('%d ',s);
            if test_mode==0
                lkt_value=lkt_ref_value_arr_as(:,mua_i,s);
            else
                lkt_value=lkt_ref_value_arr_as(:,index,s);
            end
            lkt_points_ref=lkt_value(in_place_arr);
            
            % interpolate
            if test_mode==0
                all_param_arr((mua_i-1)*size(mus_param_arr,1)+1:mua_i*size(mus_param_arr,1),2*num_layer+s)=interpn(lkt_layer_mus{1},lkt_layer_mus{2},lkt_layer_mus{4},squeeze(lkt_points_ref),mus_param_arr(:,1),mus_param_arr(:,2),mus_param_arr(:,4),'spline');
%                 all_param_arr((mua_i-1)*size(mus_param_arr,1)+1:mua_i*size(mus_param_arr,1),2*num_layer+s)=interpn(lkt_layer_mus{4},lkt_points_ref,mus_param_arr(:,4),'spline');
                all_param_arr((mua_i-1)*size(mus_param_arr,1)+1:mua_i*size(mus_param_arr,1),1:num_layer)=ones(size(mus_param_arr,1),1).*mua_param_arr(mua_i,:);
                all_param_arr((mua_i-1)*size(mus_param_arr,1)+1:mua_i*size(mus_param_arr,1),num_layer+1:2*num_layer)=mus_param_arr;
            elseif test_mode>=1
                all_param_arr(mua_i,2*num_layer+s)=interpn(lkt_layer_mus{1},lkt_layer_mus{2},lkt_layer_mus{4},squeeze(lkt_points_ref),mus_param_arr(mua_i,1),mus_param_arr(mua_i,2),mus_param_arr(mua_i,4),'spline');
%                 all_param_arr(mua_i,2*num_layer+s)=interpn(lkt_layer_mus{4},lkt_points_ref,mus_param_arr(mua_i,4),'spline');
                all_param_arr(mua_i,1:2*num_layer)=[mua_param_arr(mua_i,1:num_layer) mus_param_arr(mua_i,1:num_layer)];
            end
        end
        fprintf('\n');
    end
%     all_param_arr(:,2*num_layer+1:end)=all_param_arr(:,2*num_layer+1:end).*(setting_r./lkt_sim_set.detector_r').^2;
    save(fullfile(output_dir,'all_param_arr.mat'),'all_param_arr','-v7.3');

    if test_mode>=1
        to_save=all_param_arr(:,9:end);
        save(fullfile(output_dir,'lkt_smooth_forward.txt'),'to_save','-ascii','-tabs');
    end
end

disp('Done!');

% hist(mua_param_arr(:));

