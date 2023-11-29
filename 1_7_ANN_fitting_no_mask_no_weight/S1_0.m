%{
Find the initial value to generate database
if every subject use same parameter range of CW and TR, the result can be shared within subjects

Ting-Yi Kuo
Last update: 2023/11/1
%}

clc;clear;close all;

global lambda fitting_wl_tr Lbound Ubound cw_param_range tr_param_range

%% param
model_dir='model_arrange'; % the folder of the models

fitting_wl=(700:4:900)'; % the wavelenth to generate the spectrum
fitting_wl_tr=810;
fitting_wl_file='fitting_wl.txt'; % the file in the epsilon folder containing the fitting wavelength


% the range for the parameters
%                      hc_1    sto2_1  hc_2    sto2_2 hc_4     sto2_4  mel
mua_param_Lbound=     [1       0.3     1       0.3    50       0.3     0];
mua_param_Ubound=     [150     1       150     1      200      1       0.003 ];


%% init
wl_used_to_fitting=load(fullfile('epsilon',fitting_wl_file));
wl_to_check=[wl_used_to_fitting; fitting_wl; fitting_wl_tr];
wl_to_check=unique(wl_to_check);
lambda=wl_to_check;

ratio={[0.5],[0.33 0.66],[0.25 0.5 0.75],[0.2 0.4 0.6 0.8],[0.125 0.25 0.375 0.5 0.625 0.75 0.875]};
use_ratio=[4 3 4 3 4 3 5 3 5 3 5 3 2];


% load AK range
A_Krange=load(fullfile(model_dir,['A_Krange_arr_ZJ.mat']));
Lbound=[A_Krange.Lbound mua_param_Lbound]; % add the bound for A, K for L1, 2, 4
Ubound=[A_Krange.Ubound mua_param_Ubound];

% prepare init value
init_value=cell(1,length(Lbound));

for i=1:length(Lbound)
   ratio_select=ratio{use_ratio(i)};
   for j=1:length(ratio_select)
      temp_value=ratio_select(j)*(Ubound(i)-Lbound(i))+Lbound(i);
      init_value{i}=[init_value{i} temp_value];
   end
end

% load ANN model
load(fullfile(model_dir,['KB_cw_model.mat'])); % cw_net, cw_param_range
load(fullfile(model_dir,['ZJ_tr_model.mat'])); % tr_net, tr_param_range

%% generate the init value array
fprintf('Generating init value: \n');

fun_init_param_to_mu_spec(); % load the epsilon for target wl

init_value_arr=[];
init_value_arr_out_of_bound=[];
init_OP_arr=[];
step=1;
% total_steps=1;
% for i=1:length(use_ratio)
%     total_steps=total_steps*use_ratio(i);
% end
total_steps=4*3*4*3*4*3*7*3*7*3*7*3*2;

tic;
% progress('_start');
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
                                                    temp_init_value=[init_value{1}(i) init_value{2}(j) init_value{3}(k) init_value{4}(l) init_value{5}(m) init_value{6}(n) init_value{7}(o) init_value{8}(p) init_value{9}(q) init_value{10}(r) init_value{11}(s) init_value{12}(t) init_value{13}(u)];
                                                    [temp_OP_arr,~]=fun_param_to_mu(temp_init_value,0);
                                                    if fun_in_OP_range(temp_OP_arr)
                                                        init_value_arr=[init_value_arr; temp_init_value];
                                                        init_OP_arr(:,:,end+1)=temp_OP_arr;
                                                    else
                                                        init_value_arr_out_of_bound=[init_value_arr_out_of_bound; temp_init_value];
                                                    end
                                                    step=step+1;
%                                                     progress(step,total_steps);
                                                    if rem(step,5000)==0
                                                        fprintf('Fninsh %.2f%%: %d/%d\n',(step/total_steps)*100,step,total_steps);
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
end

process_time=toc;

save(fullfile(model_dir,'init_value_arr.mat'),'init_value_arr');

progress('_end');