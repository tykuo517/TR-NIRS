%{
Find the initial value to generate database
if every subject use same lambda and parameter range of CW and TR, the result can be shared within subjects

Ting-Yi Kuo
Last update: 2023/11/1
%}

clc;clear;close all;

global lambda fitting_wl_tr Lbound Ubound cw_param_range tr_param_range

%% param
model_dir='model_arrange'; % the folder of the models

num_init_spec=200000; % number of random generated spectrum

fitting_wl=(700:4:900)'; % the wavelenth to generate the spectrum
fitting_wl_tr=810;
fitting_wl_file='fitting_wl.txt'; % the file in the epsilon folder containing the fitting wavelength

output_dir='initValue_database_2'; % the folder to save the output database

% the range for the parameters
%                      hc_1    sto2_1  hc_2    sto2_2 hc_4     sto2_4  mel
mua_param_Lbound=     [1       0.3     1       0.3    50       0.3     0];
mua_param_Ubound=     [150     1       150     1      200      1       0.003 ];


%% init
mkdir(output_dir);

rng('shuffle'); % random seed

wl_used_to_fitting=load(fullfile('epsilon',fitting_wl_file));
wl_to_check=[wl_used_to_fitting; fitting_wl; fitting_wl_tr];
wl_to_check=unique(wl_to_check);
lambda=wl_to_check;


% load AK range
A_Krange=load(fullfile(model_dir,['A_Krange_arr_ZJ.mat']));
Lbound=[A_Krange.Lbound mua_param_Lbound]; % add the bound for A, K for L1, 2, 4
Ubound=[A_Krange.Ubound mua_param_Ubound];

% % prepare init value
% init_value=cell(1,length(Lbound));
% 
% for i=1:length(Lbound)
%    ratio_select=ratio{use_ratio(i)};
%    for j=1:length(ratio_select)
%       temp_value=ratio_select(j)*(Ubound(i)-Lbound(i))+Lbound(i);
%       init_value{i}=[init_value{i} temp_value];
%    end
% end

% load ANN model
% load(fullfile(model_dir,['ZJ_cw_model.mat'])); % cw_net, cw_param_range
% load(fullfile(model_dir,['ZJ_tr_model.mat'])); % tr_net, tr_param_range

%% generate the init value array
fprintf('Generating init value: \n');

fun_init_param_to_mu_spec(); % load the epsilon for target wl

init_value_arr=zeros(num_init_spec,size(Lbound,2));
for init_i=1:num_init_spec

    fprintf('Finding init value %d, ',init_i);

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

save(fullfile(output_dir,'init_value_arr.mat'),'init_value_arr');
