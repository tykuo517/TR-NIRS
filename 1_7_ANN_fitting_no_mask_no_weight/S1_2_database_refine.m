%{
Check if there are init value which will exceed the OP UB and LB while using the fitting wl, and then re-generate them

Benjamin Kao
Last update: 2021/01/19
%}

clc;clear;close all;clearvars -global;

global lambda fitting_wl_tr Lbound Ubound cw_net cw_param_range tr_net tr_param_range;

%% param
model_dir='model_arrange'; % the folder of the models
subject_name_arr={'KB'}; % the name of the subjects
num_SDS=7; % the number of SDS, in the ANN
num_SDS_tr=5;

num_gate=10;

num_init_spec=100000; % number of random generated spectrum

output_dir='initValue_database_2'; % the folder to save the output database

fitting_wl_file='fitting_wl.txt'; % the file in the epsilon folder containing the fitting wavelength
fitting_wl_tr=810;

% the range for the parameters
%                      hc_1    sto2_1  hc_2    sto2_2 hc_4     sto2_4  mel
mua_param_Lbound=     [1       0.3     1       0.3    50       0.3     0];
mua_param_Ubound=     [150     1       150     1      200      1       0.003 ];


%% init
rng('shuffle'); % random seed

wl_used_to_fitting=load(fullfile('epsilon',fitting_wl_file));

%% main
for sbj_i=1:length(subject_name_arr)
    
    % load ANN model
    load(fullfile(model_dir,[subject_name_arr{sbj_i} '_cw_model.mat'])); % net, param_range
    load(fullfile(model_dir,[subject_name_arr{sbj_i} '_tr_model.mat'])); % net, param_range
    
    % load AK range
    A_Krange=load(fullfile(model_dir,['A_Krange_arr_' subject_name_arr{sbj_i} '.mat']));
    Lbound=[A_Krange.Lbound mua_param_Lbound]; % add the bound for A, K for L1, 2, 4
    Ubound=[A_Krange.Ubound mua_param_Ubound];
    
    % load the already generated init DB
    init_DB=load(fullfile(output_dir,[subject_name_arr{sbj_i} '_DB.mat']));
    
    fitting_wl=init_DB.fitting_wl;
    init_spec_arr=init_DB.init_spec_arr;
    init_dtof_arr=init_DB.init_dtof_arr;
    init_value_arr=init_DB.init_value_arr;
    
    % find all the wl to check
    wl_to_check=[fitting_wl]; % [wl_used_to_fitting; fitting_wl];
    wl_to_check=unique(wl_to_check);
    
    lambda=wl_to_check;
    fun_init_param_to_mu_spec(); % load the epsilon for target wl
    
    num_re_gen=0; % number of re-generate spec
    had_retried=zeros(1,num_init_spec);
    
    for init_i=1:num_init_spec

        [OP_arr,~]=fun_param_to_mu(init_value_arr(init_i,:),0);
        
        if ~fun_in_OP_range(OP_arr)
            fprintf('\tRetry init value %d, ',init_i);
            num_re_gen=num_re_gen+1;
            had_retried(init_i)=1;
            
            %% re-generate the init value
            temp_init_value=rand(size(Lbound)).*(Ubound-Lbound)+Lbound; % generate the init param
            % adjust the K range by the A value
            for L=1:3
                temp_Krange=interp1(A_Krange.A_Krange_arr{L}(:,1),A_Krange.A_Krange_arr{L}(:,2:3),temp_init_value(2*L-1),'pchip');
                temp_init_value(2*L)=rand(1,1).*(temp_Krange(1)-temp_Krange(2))+temp_Krange(2);
            end
            [OP_arr,~]=fun_param_to_mu(temp_init_value,0);

            retry_count=1;
            while ~fun_in_OP_range(OP_arr)
                fprintf('%d ',retry_count); retry_count=retry_count+1;
                temp_init_value=rand(size(Lbound)).*(Ubound-Lbound)+Lbound;
                for L=1:3
                    temp_Krange=interp1(A_Krange.A_Krange_arr{L}(:,1),A_Krange.A_Krange_arr{L}(:,2:3),temp_init_value(2*L-1),'pchip');
                    temp_init_value(2*L)=rand(1,1).*(temp_Krange(1)-temp_Krange(2))+temp_Krange(2);
                end
                [OP_arr,~]=fun_param_to_mu(temp_init_value,0);
            end

            init_value_arr(init_i,:)=temp_init_value;
            fprintf('\n');
        end
        if rem(init_i,5000)==0
            fprintf('\tHad checked %d values...\n',init_i);
        end
    end
    
    assert(num_re_gen==length(find(had_retried>0)));

    % forward the spectrum
    lambda=fitting_wl;
    fun_init_param_to_mu_spec(); % load the epsilon for target wl
    
    % for each target spec
    for init_i=find(had_retried>0)
        fprintf('Subject %s spec %d\n',subject_name_arr{sbj_i},init_i);

        %% generate init value and spec
        [OP_arr,~]=fun_param_to_mu(init_value_arr(init_i,:),0);

        init_spec_arr(:,:,init_i)=fun_ANN_forward(OP_arr,0);
        
        OP_arr=interp1(lambda,OP_arr,fitting_wl_tr,'pchip');
        
        ann_dtof_=fun_ANN_forward(OP_arr,1);
        for i=1:num_SDS_tr
            ann_dtof(:,i)=ann_dtof_((i-1)*num_gate+1:(i-1)*num_gate+num_gate)';
        end
        init_dtof_arr(:,:,init_i)=ann_dtof;
    end

    if num_re_gen>0
        save(fullfile(output_dir,[subject_name_arr{sbj_i} '_DB2.mat']),'init_spec_arr','init_dtof_arr','init_value_arr','fitting_wl','fitting_wl_tr');
    end
    
    fprintf('Regen %d specs for %s\n',num_re_gen,subject_name_arr{sbj_i});
end

disp('Done!');