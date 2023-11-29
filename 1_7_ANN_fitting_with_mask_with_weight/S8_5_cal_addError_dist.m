%{
Calculate the dist of noise added to the fitted OPs

Benjamin Kao
Last update: 2021/03/18
%}

clc;clear;close all;

global lambda net param_range 

%% param
fitting_dir_arr={'test_fitting_2021-01-17-16-47-46','test_fitting_2021-01-22-03-42-40'};
fitting_SDS_dir_arr={'fitting_SDS1234','fitting_SDS1234','fitting_SDS346','fitting_SDS2345','fitting_SDS1234','fitting_SDS12345','fitting_SDS123456'};
fitting_dir_index=[1 1 2 2 2 2 2];
add_noise_arr=[0 1 1 1 1 1 1];
fitting_name_arr={'SDS 1234 merged no noise','SDS 1234 add system noise','SDS 346 error','SDS 2345 error','SDS 1234 error','SDS 12345 error','SDS 123456 error'};
target_name_arr={'tc_mean','ww_mean','wh_mean','yf_mean','kb_mean'}; % the name of the target spectrum
subject_fitting_SDS_arr=[7 3 5 6 5];

input_dir='fitted_result';
output_dir='saved_OPs_random';
do_plot=1;
num_random_OP=30;

to_error_OP=[1 2 3 4 7 8]; % the OPs needs to add error, 1 = mua1, 2 = mus1, 3 = mua2 ......

%% main
for sbj=1:length(target_name_arr)
    op_info=load(fullfile(input_dir,output_dir,[target_name_arr{sbj} '_OP_info.mat']));
    op_change_arr=op_info.OP_change_arr(2:1+num_random_OP,to_error_OP);
    op_change_arr=op_change_arr-1;
    
    i=subject_fitting_SDS_arr(sbj);
    if add_noise_arr(i)==0
        to_load_name='OP_error_arr_noError.mat';
    else
        to_load_name='OP_error_arr_Error.mat';
    end
    error_info=load(fullfile(fitting_dir_arr{fitting_dir_index(i)},'arrangement',fitting_SDS_dir_arr{i},to_load_name));
    error_pdf_arr=error_info.OP_pdf_arr;
    error_pdf_binCenter=error_info.OP_pdf_binCenter;
    
    OP_name_arr={'\mu_{a,scalp}','\mu_{s,scalp}','\mu_{a,skull}','\mu_{s,skull}','\mu_{a,GM}','\mu_{s,GM}'};
%     figure('Units','pixels','position',[0 0 1920 1080]);
    figure('Units','pixels','position',[0 0 1400 720]);
    ti=tiledlayout('flow','TileSpacing','compact','Padding','none');
    for opi=1:6
        [a,b]=hist(op_change_arr(:,opi),50);
%         a=a./sum(a);

        nexttile();
%         plot(error_pdf_binCenter{opi}*100,error_pdf_arr{opi},b*100,a,'-.','LineWidth',2);
        hold on;
        plot(error_pdf_binCenter{opi}*100,error_pdf_arr{opi},'LineWidth',2);        
        ylabel('PDF');
        yyaxis right
        bar(b*100,a,1);
        ylabel('sampled error distribution');
        xlabel([OP_name_arr{opi} ' error (%)']);
        grid on;
        legend({'original error PDF','sampled error dist'},'Location','northeast');
        set(gca,'fontsize',12, 'FontName', 'Times New Roman');
    end
    print(fullfile(input_dir,output_dir,[target_name_arr{sbj} '_OP_dist.png']),'-dpng','-r300');
    close all;
end

disp('Done!');