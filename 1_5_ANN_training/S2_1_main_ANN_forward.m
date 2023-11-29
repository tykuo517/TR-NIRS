%{
Use the trained ANN to generate spec, in order to compare with the simulated spectrum or lut spectrum

Benjamin Kao
Last update: 2020/12/23
%}

clc;clear;close all;

%% param
% model_dir='ZJ_2020-12-22-12-36-51'; % the model train folder
model_dir_arr=load('model_dir.mat'); % the model train folder

for i=1:size(model_dir_arr.model_dir,1)
    clearvars -except i model_dir_arr
    
    model_dir=model_dir_arr.model_dir{i,2};
    
    %% test
    fun_ANN_init(model_dir);
    
    %% test 1
    op=load('OPs_to_sim_6/toSim_OP_1.txt');
    op=op(:,1:8);
    spec=fun_ANN_forward(op);
    save(fullfile(model_dir,['ANN_forward_test1.txt']),'spec','-ascii','-tabs');
    
    %% test 2
    op=load('OPs_to_sim_11/toSim_OP_65.txt');
    op=op(:,1:8);
    spec=fun_ANN_forward(op);
    save(fullfile(model_dir,['ANN_forward_test2.txt']),'spec','-ascii','-tabs');
    
    %% test 3
    op=load('OPs_to_sim_11/toSim_OP_66.txt');
    op=op(:,1:8);
    spec=fun_ANN_forward(op);
    save(fullfile(model_dir,['ANN_forward_test3.txt']),'spec','-ascii','-tabs');
end

disp('Done!');