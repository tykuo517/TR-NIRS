%{
Save the error PDF for each subject

Benjamin Kao
Last update: 2021/03/17
%}

clc;clear;close all; clearvars -global;

%% param
output_dir='fitted_result';
fitting_dir_arr={'test_fitting_2021-01-17-16-47-46','test_fitting_2021-01-22-03-42-40'}; % the mother dir of the testing target fitting
fitting_SDS_dir_arr={'fitting_SDS1234','fitting_SDS1234','fitting_SDS346','fitting_SDS2345','fitting_SDS1234','fitting_SDS12345','fitting_SDS123456'}; % the sub fitting dir in the testing target fitting dir for each SDS combination
fitting_dir_index=[1 1 2 2 2 2 2];
add_noise_arr=[0 1 1 1 1 1 1];
fitting_name_arr={'SDS 1234 merged no noise','SDS 1234 add system noise','SDS 346 error','SDS 2345 error','SDS 1234 error','SDS 12345 error','SDS 123456 error'}; % the name for each SDS or error condition, just for understanding
target_name_arr={'tc_mean','ww_mean','wh_mean','yf_mean','kb_mean'}; % the name of the target spectrum
subject_fitting_SDS_arr=[7 3 5 6 5];
IS_n_interval=10000; % how many inverse sampling point on the PDF function
do_test_IS=1;

%% main
for sbj=1:length(target_name_arr)
    i=subject_fitting_SDS_arr(sbj);
    if add_noise_arr(i)==0
        to_load_name='OP_error_arr_noError.mat';
    else
        to_load_name='OP_error_arr_Error.mat';
    end
    error_info=load(fullfile(fitting_dir_arr{fitting_dir_index(i)},'arrangement',fitting_SDS_dir_arr{i},to_load_name));
    error_pdf_arr=error_info.OP_pdf_arr;
    error_pdf_binCenter=error_info.OP_pdf_binCenter;
    
    error_CDF_arr={};
    error_IS_arr={}; % inverse transform sampling
    for opi=1:6
        error_CDF_arr{opi}=cumsum(error_pdf_arr{opi});
        [unique_CDF,unique_index]=unique(error_CDF_arr{opi});
        CDF_points=error_pdf_binCenter{opi};
        CDF_points(2:end)=CDF_points(2:end)+(CDF_points(3)-CDF_points(2))/2;
        CDF_points(1)=CDF_points(1)-(CDF_points(3)-CDF_points(2))/2;
        CDF_points=CDF_points(unique_index);
        error_IS_arr{opi}=interp1(unique_CDF,CDF_points,linspace(0,1,IS_n_interval+1),'pchip');
    end
    
    % test IS
    if do_test_IS
        OP_name_arr={'\mu_{a,scalp}','\mu_{s,scalp}','\mu_{a,skull}','\mu_{s,skull}','\mu_{a,GM}','\mu_{s,GM}'};
%         figure('Units','pixels','position',[0 0 1920 1080]);
        figure('Units','pixels','position',[0 0 1280 720]);
        ti=tiledlayout('flow','TileSpacing','compact','Padding','none');
        n_test=100000;
        for opi=1:6
            test_sample_point=interp1(linspace(0,1,IS_n_interval+1),error_IS_arr{opi},rand(n_test,1),'pchip');
            [a,b]=hist(test_sample_point,100);
            a=a./sum(a);

            nexttile();
            plot(error_pdf_binCenter{opi}*100,error_pdf_arr{opi},b*100,a,'-.','LineWidth',2);
            xlabel([OP_name_arr{opi} ' error (%)']);
            ylabel('PDF');
            legend({'original error PDF','sampled error PDF'},'Location','best');
            set(gca,'fontsize',12, 'FontName', 'Times New Roman');
        end
        print(fullfile(output_dir,[target_name_arr{sbj} '_error_IS.png']),'-dpng','-r300');
        close all;
    end
    
    save(fullfile(output_dir,[target_name_arr{sbj} '_error_IS_arr.mat']),'error_IS_arr');
end

disp('Done!');