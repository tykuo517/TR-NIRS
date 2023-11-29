%{
Check if there are init value which will exceed the OP UB and LB while using the fitting wl, and then re-generate them

Benjamin Kao
Last update: 2021/01/19
%}

clc;clear;close all;clearvars -global;

%% param
subject_name_arr={'ZJ','WW2','YF','YH','WH2','KB','SJ','BT','SC'}; % the name of the subjects
num_SDS=7; % the number of SDS, in the ANN

num_init_spec=50000; % number of random generated spectrum

output_dir='initValue_database_2'; % the folder to save the output database

for sbj=1:length(subject_name_arr)
    orig_DB=load(fullfile(output_dir,[subject_name_arr{sbj} '_DB.mat']));
    new_DB=load(fullfile(output_dir,[subject_name_arr{sbj} '_DB2.mat']));
    
    init_value_diff=sum(orig_DB.init_value_arr~=new_DB.init_value_arr,2)>0;
    init_value_diff_index=find(init_value_diff);
    
    init_spec_diff=sum(orig_DB.init_value_arr~=new_DB.init_value_arr,2:3)>0;
    init_spec_diff_index=find(init_spec_diff);
    
    assert(length(init_value_diff_index)==length(init_spec_diff_index));
    assert(sum(init_value_diff_index~=init_spec_diff_index)==0);
    fprintf('There are %d diffenent init value in %s.\n',length(init_spec_diff_index),subject_name_arr{sbj});
    DoDelete=input('Replace the old file? 1 or 0?');
    if DoDelete==1
        delete(fullfile(output_dir,[subject_name_arr{sbj} '_DB.mat']));
        movefile(fullfile(output_dir,[subject_name_arr{sbj} '_DB2.mat']),fullfile(output_dir,[subject_name_arr{sbj} '_DB.mat']));
    end
end

disp('Done!');