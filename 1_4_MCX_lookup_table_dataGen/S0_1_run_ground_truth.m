 %{
Run 2E11 photons as ground truth for smoothing analysis.
Change the param of scalp, skull, CSF and GM, mus of WM = 3* GM

Ting-Yi Kuo
Last update: 2023/10/16
Version: 4.41
%}

clc;clear;close all;
%% param
model_folder='models'; % the folder of the models
target_name_arr={'KB'}; % the name of the model to simulate
num_photon=2E11; % the number of photon
num_SDS=5; % number of detectors
layer_mus={[100,225],[50,200],23,[25:25:275]}; % mus for each layer, 1/cm
n=1.4;
g=0.9;

sim_version=4.41;

%% init
% make the mus table
mus_table=[];
for i=1:length(layer_mus{1})
    for j=1:length(layer_mus{3})
        for k=1:length(layer_mus{4})
            mus_table=[mus_table; layer_mus{1}(i) layer_mus{2}(i) layer_mus{3}(j) layer_mus{4}(k) layer_mus{4}(k)*3];
        end
    end
end

    
% simulation setting
sim_set.num_SDS=num_SDS;
sim_set.detector_r=[1.5 1.5 1.5 1.5 1.5]'; % mm
sim_set.detector_larger_r=sim_set.detector_r.*2;
sim_set.detector_NA=ones(sim_set.num_SDS,1)*0.26;
sim_set.fiber_n=1.457;
sim_set.num_photon=num_photon;
sim_set.photon_per_simulation=1E8;%200000000;
sim_set.mcx_max_detpt=6E7;
sim_set.source_type='cone';
sim_set.source_NA=0.26;
sim_set.source_r=0;


%% main
for ti=1:length(target_name_arr)
    clearvars -except target_name_arr ti model_folder layer_mus mus_table sim_set sim_version n g;
    
    fprintf('Processing subject %d/%d\n',ti,length(target_name_arr));
    
    %% init
    target_name=target_name_arr{ti};
    
    target_folder=target_name;
    mkdir(fullfile(target_folder));
    
    fid=fopen(fullfile(target_folder,'sim_version.txt'),'w');
    fprintf(fid,'%.2f',sim_version);
    fclose(fid);
    
    MRI_voxel_file=fullfile(model_folder,['headModel' target_name '_EEG.mat']);
    load(MRI_voxel_file);
    
    if exist('model_version','var')==0
        model_version=1;
        do_sinus=0;
    end
    
    % change some setting according to model
    if do_sinus==1
        sim_set.num_layer=6;
        sim_set.to_output_layer=1:6;
        sim_set.n=[n n n n n 1];
        sim_set.g=[g g g g g 1];
    else
        sim_set.num_layer=5;
        sim_set.to_output_layer=1:5;
        sim_set.n=ones(sim_set.num_layer,1)*n;
        sim_set.g=ones(sim_set.num_layer,1)*g;
    end
    
    %% save the setting
    save(fullfile(target_folder,'layer_mus.mat'),'layer_mus');
    
    %% save the settings
    copyfile(fullfile(model_folder,[target_name '_probe_pos.txt']),fullfile(target_folder,[target_name '_probe_pos.txt']));
    copyfile(fullfile(model_folder,[target_name '_probe_dir.txt']),fullfile(target_folder,[target_name '_probe_dir.txt']));
    save(fullfile(target_folder,'mus_table.txt'),'mus_table','-ascii','-tabs');
    save(fullfile(target_folder,'sim_set.mat'),'sim_set');
    
    % backup this script
    copyfile([mfilename '.m'],fullfile(target_folder,[mfilename '.m']));
    
    %% final
    fprintf('There will be %d simulation.\n',size(mus_table,1));
    
    % do simulation
    fun_MCX_run_lookup(target_name);
    
    fprintf('Done!');
end