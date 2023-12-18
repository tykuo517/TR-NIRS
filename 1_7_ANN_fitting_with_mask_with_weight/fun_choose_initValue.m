%{
Using the database to choose the initial value of the fitting

Benjamin Kao
Last update: 2021/01/12
%}

function [sorted_param_arr,sorted_cw_error,sorted_tr_error,sorted_error]=fun_choose_initValue(target_spec,target_wl,target_dtof,target_name)
global cw_flag SDS_choosed_cw SDS_sim_correspond_cw SDS_choosed_tr tr_flag SDS_sim_correspond_tr gate_sim_correspond mask weight
%% param
dataBase_dir='initValue_database_2';

%% load
init_DB=load(fullfile(dataBase_dir,[target_name '_DB.mat']));
% CW
if cw_flag
    target_spec=interp1(target_wl,target_spec,init_DB.fitting_wl);
    error_arr=init_DB.init_spec_arr(:,SDS_sim_correspond_cw,:)./target_spec-1;
    error_arr=error_arr(:,SDS_choosed_cw,:);
    total_error_arr(:,1)=squeeze(sqrt(mean(mean(error_arr.^2,1),2)));
    total_error_arr(:,1)=total_error_arr(:,1).*weight(1)./sum(weight);
else
    total_error_arr(:,1)=NaN*ones(1,size(init_DB.init_spec_arr,3));
end

% TR
if tr_flag
    error_arr=init_DB.init_dtof_arr(gate_sim_correspond,SDS_sim_correspond_tr,:)./target_dtof-1;
    error_arr=error_arr.*mask; % choose gate
    error_arr=error_arr(:,SDS_choosed_tr,:); % choose SDS
    select_gate=sum(mask);
    total_error_arr(:,2)=squeeze(sqrt(mean(sum(error_arr.^2,1)./select_gate(SDS_choosed_tr),2)));
    total_error_arr(:,2)=total_error_arr(:,2).*weight(2)./sum(weight);
else
    total_error_arr(:,2)=NaN*ones(1,size(init_DB.init_dtof_arr,3));
end

% calculate total error
if sum(isnan(total_error_arr))==0
    total_error=sum(total_error_arr,2);
else
    total_error=sum(total_error_arr(~isnan(total_error_arr)),2);
end
[sorted_error, sorted_index]=sort(total_error);

cw_error=total_error_arr(:,1);
sorted_cw_error=cw_error(sorted_index,:);
tr_error=total_error_arr(:,2);
sorted_tr_error=tr_error(sorted_index,:);

sorted_param_arr=init_DB.init_value_arr(sorted_index,:);
end



