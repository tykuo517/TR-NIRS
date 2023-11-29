

clear;

input_dir='test_fitting_2023-11-02-11-10-06'; % please move the fitting folders into this folder first.
subject_name_arr={'KB'};
num_anser_to_generate=15; % number of target spec (true answer)
num_error_to_generate=15; % number of adding noise to the same, the first one will have no error
times_to_fitting=20; % number of fitting using different init value

fitting_dir='fitting_SDS_234_gw'; % the fitting folder

for sbj=1:length(subject_name_arr)

    for target_i=1:num_anser_to_generate
        for error_i=1:num_error_to_generate

            temp_target_folder=fullfile(input_dir,subject_name_arr{sbj},['target_' num2str(target_i) '_' num2str(error_i)],fitting_dir);
            
%             if exist(fullfile(temp_target_folder,'fitting_results.csv'),'file')
%                 delete(fullfile(temp_target_folder,'fitting_results.csv'));
%             end

%             if exist(temp_target_folder,'file')
%                 rmdir(temp_target_folder,'s');
%             end

            newFolderName = fullfile(input_dir,subject_name_arr{sbj},['target_' num2str(target_i) '_' num2str(error_i)],[fitting_dir '_3']);  % 你的新文件夹名
            status = movefile(temp_target_folder, newFolderName);
            
        end
    end
end