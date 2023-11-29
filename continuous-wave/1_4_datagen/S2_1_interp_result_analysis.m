%{
Plot reflectance before amd after interpolation

Ting-Yi Kuo 
Last update: 2023/10/19
%}

clc;clear;close all;

%% param
subject_name_arr={'WH'}; % the name of the subject
input_dir_arr={'WH_2023-11-14-21-53-01'}; % the lut folder of each subject
interp_dir_arr={'WH_test1_2023-11-14-21-57-04'};

for sbj_i=1:length(subject_name_arr)
    
    subject_name=subject_name_arr{sbj_i};
    input_dir=input_dir_arr{sbj_i};
    interp_dir=interp_dir_arr{sbj_i};
    num_SDS=5;
    num_gate=10;
    
    % lookup table value
    all_param_arr_lookup=load(fullfile(input_dir,'all_param_arr.mat'));
    all_param_arr_lookup=all_param_arr_lookup.all_param_arr;
    index=all(all_param_arr_lookup(:,1:7)==[0.45 0.3 0.042 0.4 100 50 23], 2); % set the mua,mus combination you want to plot 
    param_arr=all_param_arr_lookup(index,:);
    
    % interpolation value
    param_arr_interp=load(fullfile(interp_dir,'all_param_arr.mat'));
    param_arr_interp=param_arr_interp.all_param_arr;
    
    figure('Units','pixels','position',[0 0 1920 1080]);
    ti=tiledlayout(4,6);
    for s=2:num_SDS
        for g=1:6 %num_gate
            nexttile;
            plot(param_arr(:,8),param_arr(:,8+g+(s-1)*num_gate),'-o');
            hold on
            plot(param_arr_interp(:,8),param_arr_interp(:,8+g+(s-1)*num_gate));
            title(['SDS ' num2str(s) ' Gate ' num2str(g)]);
            xlabel('\mu_{s,GM}');
            ylabel('reflectance');
        end
    end
    title(ti,subject_name);
    print(fullfile('results_interp',['interp_b_and_a_' subject_name '.png']),'-dpng','-r200');
    
end

