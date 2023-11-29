%{
Load and init the ANN model

Benjamin Kao
Last update: 2020/10/26
%}

function fun_ANN_init(input_dir)
global net param_range;
fprintf('Loading ANN from %s\n',input_dir);
load(fullfile(input_dir,'ANN_model.mat')); % net

param_range=load(fullfile(input_dir,'param_range.txt'));

end