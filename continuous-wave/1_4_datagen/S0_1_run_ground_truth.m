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
lookup_table_arr='../1_3_MCX_lookup_table_CW'; % the dir containing the unmerged lookup table
target_name_arr={'KB'}; % the name of the model to simulate
num_photon=1E11; % the number of photon
n=1.4;
g=0.9;

sim_version=4.41;


%% main
for ti=1:length(target_name_arr)
    
    fprintf('Processing subject %d/%d\n',ti,length(target_name_arr));
    
    %% init
    target_name=target_name_arr{ti};
    lkt_dir=fullfile(lookup_table_arr,target_name); % the dir containing the unmerged lookup table
    
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
    
    % load lookup table information
    sim_set=load(fullfile(lkt_dir,'sim_set.mat'));
    sim_set=sim_set.sim_set;
    lkt_layer_mus=load(fullfile(lkt_dir,'layer_mus.mat'));
    lkt_layer_mus=lkt_layer_mus.layer_mus;
    lkt_mus_table=load(fullfile(lkt_dir,'mus_table.txt'));
         
    sim_set.num_photon=num_photon;
    sim_set.photon_per_simulation=1E9;
    sim_set.mcx_max_detpt=6E7;
    
    
    % make the mus table
    layer_mus={lkt_layer_mus{1}([1 7 end]),lkt_layer_mus{2}([1 5 end]),lkt_layer_mus{3}([1 2 end]),lkt_layer_mus{4}}; % mus for each layer, 1/cm
    mus_table=[];
    for i=1:length(layer_mus{1})    
        for j=1:length(layer_mus{4})
            mus_table=[mus_table; layer_mus{1}(i) layer_mus{2}(i) layer_mus{3}(i) layer_mus{4}(j) layer_mus{4}(j)*3];
        end
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