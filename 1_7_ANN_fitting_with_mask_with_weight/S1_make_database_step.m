%{
Generate the database for the fitting program to find the initial value

Benjamin Kao
Last update: 2021/01/19
%}

clc;clear;close all;clearvars -global;

global lambda fitting_wl_tr Lbound Ubound cw_net cw_param_range tr_net tr_param_range;

%% param
model_dir='model_arrange'; % the folder of the models
subject_name_arr={'KB'}; % the name of the subjects
num_SDS_cw=7;
num_SDS_tr=5;
num_gate=10; % the number of SDS, in the ANN

num_init_spec=150000; % number of random generated spectrum
ratio={[0.5],[0.33 0.66],[0.25 0.5 0.75],[0.2 0.4 0.6 0.8]};
use_ratio=[3 2 3 2 3 2 4 2 4 2 4 2 2];


fitting_wl=(700:4:900)'; % the wavelenth to generate the spectrum
fitting_wl_tr=810;

fitting_wl_file='fitting_wl.txt'; % the file in the epsilon folder containing the fitting wavelength

output_dir='initValue_database_2'; % the folder to save the output database

% the range for the parameters
%                      hc_1    sto2_1  hc_2    sto2_2 hc_4     sto2_4  mel
mua_param_Lbound=     [1       0.3     1       0.3    50       0.3     0];
mua_param_Ubound=     [150     1       150     1      200      1       0.003 ];

% %                      hc_2    sto2_2  hc_4  sto2_4 
% mua_param_Lbound=     [1       0.3     50      0.3];
% mua_param_Ubound=     [150     1       200     1 ];


%% init
mkdir(output_dir);

rng('shuffle'); % random seed

wl_used_to_fitting=load(fullfile('epsilon',fitting_wl_file));

% find all the wl to check
wl_to_check=[wl_used_to_fitting; fitting_wl];
wl_to_check=unique(wl_to_check);

%% main

for sbj_i=1:length(subject_name_arr)
    lambda=wl_to_check;
    fun_init_param_to_mu_spec(); % load the epsilon for target wl
    
    % load ANN model
    load(fullfile(model_dir,[subject_name_arr{sbj_i} '_cw_model.mat'])); % cw_net, cw_param_range
    load(fullfile(model_dir,[subject_name_arr{sbj_i} '_tr_model.mat'])); % tr_net, tr_param_range
    
    % load AK range
    A_Krange=load(fullfile(model_dir,['A_Krange_arr_' subject_name_arr{sbj_i} '.mat']));
    Lbound=[A_Krange.Lbound mua_param_Lbound]; % add the bound for A, K for L1, 2, 4
    Ubound=[A_Krange.Ubound mua_param_Ubound];


    % generate the init value array
%     init_value_arr=zeros(num_init_spec,size(Lbound,2));

        
    %% generate init value
    init_value=cell(1,length(Lbound));

    for i=1:length(Lbound)
       ratio_select=ratio{use_ratio(i)};
       for j=1:length(ratio_select)
          temp_value=ratio_select(j)*(Ubound(i)-Lbound(i))+Lbound(i);
          init_value{i}=[init_value{i} temp_value];
       end
    end

    fprintf('Generating init value: ');

    init_value_arr=[];
    init_value_arr_out_of_bound=[];
    OP_arr=[];
    step=1;
    total_steps=1;
    for i=1:length(use_ratio)
        total_steps=total_steps*use_ratio(i);
    end
%     totalsteps=length(ratio)^13;
    progress('_start');
%     h=waitbar(0, 'Generating init value...');
    for i=1:length(init_value{1})
        for j=1:length(init_value{2})
            for k=1:length(init_value{3})
                for l=1:length(init_value{4})
                    for m=1:length(init_value{5})
                        for n=1:length(init_value{6})
                            for o=1:length(init_value{7})
                                for p=1:length(init_value{8})
                                    for q=1:length(init_value{9})
                                        for r=1:length(init_value{10})
                                            for s=1:length(init_value{11})
                                                for t=1:length(init_value{12})
                                                    for u=1:length(init_value{13})
%                                                         fprintf('%d,',step);
                                                        temp_init_value=[init_value{1}(i) init_value{2}(j) init_value{3}(k) init_value{4}(l) init_value{5}(m) init_value{6}(n) init_value{7}(o) init_value{8}(p) init_value{9}(q) init_value{10}(r) init_value{11}(s) init_value{12}(t) init_value{13}(u)];
                                                        [temp_OP_arr,~]=fun_param_to_mu(temp_init_value,0);
                                                        if fun_in_OP_range(temp_OP_arr)
                                                            init_value_arr=[init_value_arr; temp_init_value];
                                                            init_OP_arr=[init_OP_arr;temp_OP_arr];
                                                        else
                                                            init_value_arr_out_of_bound=[init_value_arr_out_of_bound; temp_init_value];
                                                        end
                                                        step=step+1;
                                                        progress(step,total_steps);
%                                                         waitbar(step/totalsteps,h,sprintf('Generating init value...(%d/%d)',step,totalsteps));
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    progress('_end');
  

    % forward the spectrum
    lambda=fitting_wl;
    fun_init_param_to_mu_spec(); % load the epsilon for target wl
    
    init_spec_arr=zeros(length(fitting_wl),num_SDS_cw,length(init_value_arr));
    init_dtof_arr=zeros(num_gate,num_SDS_tr,length(init_value_arr));

    % for each target spec
    for init_i=1:size(init_value_arr,1)
        fprintf('Subject %d spec %d\n',sbj_i,init_i);

        %% generate init value and spec
        [OP_arr,~]=fun_param_to_mu(init_value_arr(init_i,:),0);

        init_spec_arr(:,:,init_i)=fun_ANN_forward(OP_arr,0);
        
        temp_OP_arr=[fitting_wl OP_arr];
        OP_arr=interp1(temp_OP_arr(:,1),temp_OP_arr(:,2:end),fitting_wl_tr,'pchip');
        
        target_dtof_=fun_ANN_forward(OP_arr,1);
        
        for i=1:num_SDS_tr
            target_dtof(:,i)=target_dtof_((i-1)*num_gate+1:i*num_gate)';
        end
        
        init_dtof_arr(:,:,init_i)=target_dtof;
        
    end

    save(fullfile(output_dir,[subject_name_arr{sbj_i} '_DB.mat']),'init_spec_arr','init_dtof_arr','init_value_arr','fitting_wl','fitting_wl_tr');
end

disp('Done!');