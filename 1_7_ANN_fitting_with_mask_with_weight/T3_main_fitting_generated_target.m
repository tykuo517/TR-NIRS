%{
Fitting the generated target spectrum

Benjamin Kao
Last update: 2021/01/21
%}

clc;clear;close all;

global lambda fitting_wl_tr Lbound Ubound cw_flag cw_net cw_param_range tr_flag tr_net tr_param_range mask weight;
global target_spec orig_target_spec target_dtof;
global num_SDS_cw num_SDS_tr SDS_choosed_cw SDS_choosed_tr A_Krange SDS_sim_correspond_cw SDS_sim_correspond_tr num_gate gate_choosed gate_sim_correspond;

%% param
subject_name_arr={'KB','WH','ZJ'}; %,'WH','ZJ'
num_anser_to_generate=10; % number of target spec (true answer)
num_error_to_generate=1; % number of adding noise to the same, the first one will have no error
times_to_fitting=20; % number of fitting using different init value
input_dir='test_fitting_2023-12-12-11-59-43'; % the folder containing the generated target spectrum


fitting_wl_file='fitting_wl.txt'; % the file in the epsilon folder containing the fitting wavelength
fitting_wl_tr=810;

model_dir='model_arrange';

%% CW setting
cw_flag=0;
num_SDS_cw=6; % how many SDS are in the target spectrum

if cw_flag
    SDS_choosed_cw=[1]; % the SDS chosen to fitted in the target spectrum
    SDS_sim_correspond_cw=[1 2 3 4 5 6]; % the SDS index in the simulated spectrum corresponding to each SDS in the target spectrum, SDS_sim_correspond(x)=y means 'measure x = sim y'
end

%% TR setting
tr_flag=1;
num_SDS_tr=5;
num_gate=10;

% SDS
if tr_flag
    SDS_choosed_tr=[1]; % the SDS chosen to fitted in the target spectrum 
    SDS_sim_correspond_tr=[1 2 3 4 5]; % the SDS index in the simulated spectrum corresponding to each SDS in the target spectrum, SDS_sim_correspond(x)=y means 'measure x = sim y'

    % choose gate or not
    choose_gate_flag=1; % =0: fitted the gate within 50%~0.01% max value
    num_gate=10;
    gate_choosed=[1 2]; % if choose_gate_manually==1
    gate_sim_correspond=[1 2 3 4 5 6 7 8 9 10];

    % use gate weight or not
    gate_weight_flag=0; % weight between gate
    %if use_gate_weight==1
    weight_SDS=[2 3 4];
    weight_gate=[7];
    weight_value=2;
end

%% CW and TR weight setting
sys_weight_flag=0; % weight between CW and TR


%% init
% Named output directory
output_dir='fitting';
if cw_flag
    output_dir=[output_dir '_cw'];
    for s=1:length(SDS_choosed_cw)
        output_dir=[output_dir num2str(SDS_choosed_cw(s))];
    end
end

if tr_flag
    output_dir=[output_dir '_tr'];
    for s=1:length(SDS_choosed_tr)
        output_dir=[output_dir num2str(SDS_choosed_tr(s))];
    end
    
    if choose_gate_flag
        output_dir=[output_dir '_gate'];
        for g=1:length(gate_choosed)
            output_dir=[output_dir num2str(gate_choosed(g))];
        end
    end
    
    if gate_weight_flag==1
        for w=1:length(weight_gate)
            output_dir=[output_dir '_gw_' num2str(weight_gate(w))];
        end
    end
end

if cw_flag && tr_flag && sys_weight_flag
    output_dir=[output_dir '_sw']; % weight between different system (CW and TR)
end



%                 hc_1    sto2_1 hc_2    sto2_2 hc_4     sto2_4   mel
mua_param_Lbound=[1       0.3     1       0.3    50       0.3     0];
mua_param_Ubound=[150     1       150     1      200      1       0.003 ];

% %                      hc_2    sto2_2  hc_4  sto2_4 
% mua_param_Lbound=     [1       0.3     50      0.3];
% mua_param_Ubound=     [150     1       200     1 ];

mua_coef_bound=[mua_param_Lbound; mua_param_Ubound];


fitting_wl=load(fullfile('epsilon',fitting_wl_file));
fitting_wl=[fitting_wl;fitting_wl_tr];
lambda=unique(fitting_wl);

fun_init_param_to_mu_spec(); % load the epsilon for fitting wl

% backup this script
copyfile([mfilename '.m'],fullfile(input_dir,[mfilename '.m']));


%% main
for sbj=1:length(subject_name_arr)
    % load ANN model
    load(fullfile(model_dir,[subject_name_arr{sbj} '_cw_model.mat'])); % net, param_range
    load(fullfile(model_dir,[subject_name_arr{sbj} '_tr_model.mat'])); % net, param_range
    
    % load the param bound
    A_Krange=load(fullfile(model_dir,['A_Krange_arr_' subject_name_arr{sbj} '.mat']));
    Lbound=[A_Krange.Lbound mua_param_Lbound]; % add the bound for A, K for L1, 2, 4
    Ubound=[A_Krange.Ubound mua_param_Ubound];
    
    % for each target spec
    for target_i=1:num_anser_to_generate
        % for each error spec
        for error_i=1:num_error_to_generate
            thisSpec_outputDir=fullfile(input_dir,subject_name_arr{sbj},['target_' num2str(target_i) '_' num2str(error_i)],output_dir);
            if exist(thisSpec_outputDir,'dir')==0
                mkdir(thisSpec_outputDir);
            end
            
            % check if all fitting had been done
            target_had_been_done=1;
            for fitting_i=1:times_to_fitting
                if exist(fullfile(thisSpec_outputDir,['fitting_' num2str(fitting_i)],'fitted_spec.txt'),'file')==0
                    target_had_been_done=0;
                    break;
                end
            end
            if target_had_been_done==1
                fprintf('Skip %s\n',thisSpec_outputDir);
                continue;
            end

            % load target spec, and interpolation
            orig_target_spec=load(fullfile(input_dir,subject_name_arr{sbj},['target_' num2str(target_i) '_' num2str(error_i)],'target_spec.txt'));
            target_spec=interp1(orig_target_spec(:,1),orig_target_spec(:,2:end),lambda);
            
            orig_target_dtof=load(fullfile(input_dir,subject_name_arr{sbj},['target_' num2str(target_i) '_' num2str(error_i)],'target_dtof.txt'));
            target_dtof=orig_target_dtof;
            
            % choose gate to fit
            if tr_flag
                mask=zeros(num_gate,num_SDS_tr);
                if choose_gate_flag
                    mask(gate_choosed,:)=1;
                else
                    for s=1:num_SDS_tr
                        dtof=target_dtof(:,s);
                        [max_value,index]=max(dtof);
                        for g=1:num_gate
                            if g<index
                                mask(g,s)=dtof(g)>=0.5*max_value;
                            else
                                mask(g,s)=dtof(g)>=0.0001*max_value;
                            end
                        end
                    end
                end
            end
                
            % assign weight
            if cw_flag && tr_flag && sys_weight_flag==1
                weight=[length(lambda)*length(SDS_choosed_cw) sum(mask(:,SDS_choosed_tr),'all')];
            elseif cw_flag && tr_flag && sys_weight_flag==0
                weight=[1 1];
            elseif ~cw_flag && tr_flag
                weight=[0 1];
            elseif cw_flag && ~tr_flag
                weight=[1 0];
            end
            
            % assign weight to gate (if needed edit here)
            if tr_flag && gate_weight_flag
                for s=weight_SDS
                    for g=weight_gate
                        mask(g,s)=weight_value;
                    end
                end
            end
            
            %% fitting
            % find the init value
            [init_arr,init_cw_error,init_tr_error,init_error]=fun_choose_initValue(orig_target_spec(:,2:end),orig_target_spec(:,1),target_dtof,subject_name_arr{sbj});
            init_arr=init_arr(1:times_to_fitting,:);
            init_cw_error=init_cw_error(1:times_to_fitting,:);
            init_tr_error=init_tr_error(1:times_to_fitting,:);
            
            init_error=init_error(1:times_to_fitting);
            to_save=[init_arr init_cw_error init_tr_error init_error];

            fid=fopen(fullfile(thisSpec_outputDir,'init_arr.csv'),'a');
            fprintf(fid,'A1,K1,A2,K2,A4,K4,hc_1,StO2_1,hc_2,StO2_2,hc_4,StO2_4,mel,cw_error,tr_error,error,\n');
            for fitting_i=1:times_to_fitting
                fprintf(fid,'%f,',to_save(fitting_i,:));
                fprintf(fid,'\n');
            end
            fclose(fid);
            save(fullfile(thisSpec_outputDir,'init_arr.txt'),'to_save','-ascii','-tabs');
            if cw_flag && tr_flag
                save(fullfile(thisSpec_outputDir,'fitting_info.mat'),'cw_flag','tr_flag','SDS_choosed_cw','SDS_choosed_tr','mask','weight');
            elseif cw_flag && ~tr_flag
                save(fullfile(thisSpec_outputDir,'fitting_info.mat'),'cw_flag','tr_flag','SDS_choosed_cw','mask','weight');
            elseif ~cw_flag && tr_flag
                save(fullfile(thisSpec_outputDir,'fitting_info.mat'),'cw_flag','tr_flag','SDS_choosed_tr','mask','weight');
            else
                error('Please choose at least one system to fit.');
            end

            % start fitting using different init value
            for fitting_i=1:times_to_fitting
                had_fitted=exist(fullfile(thisSpec_outputDir,['fitting_' num2str(fitting_i)],'fitted_spec.txt'),'file') | exist(fullfile(thisSpec_outputDir,['fitting_' num2str(fitting_i)],'fitted_dtof.txt'),'file');
                if had_fitted~=0
                    fprintf('Skip %s\n',fullfile(thisSpec_outputDir,['fitting_' num2str(fitting_i)]));
                end
                while had_fitted==0
                    mkdir(fullfile(thisSpec_outputDir,['fitting_' num2str(fitting_i)]));

                    fun_auto_fittingMuaMus(init_arr(fitting_i,:),fullfile(thisSpec_outputDir,['fitting_' num2str(fitting_i)]));
                    
                    had_fitted=exist(fullfile(thisSpec_outputDir,['fitting_' num2str(fitting_i)],'fitted_spec.txt'),'file') | exist(fullfile(thisSpec_outputDir,['fitting_' num2str(fitting_i)],'fitted_dtof.txt'),'file');
                end
           end
        end
    end
end

%% function