%{
Check if the OPs are in the range

Output:
result: =1 if the OP is in range, =0 if OP is out of range

Benjamin Kao
Last update: 2020/09/19
%}

function result=fun_in_OP_range(OP_arr)
global cw_param_range tr_param_range lambda fitting_wl_tr;
OP_arr=OP_arr(:,[1:2:8 2:2:8]); % mua and mus
OP_arr_tr=interp1(lambda,OP_arr,fitting_wl_tr,'spline');

result=OP_arr<cw_param_range(2,:) | OP_arr>cw_param_range(1,:); % if the OP is out of range
temp_result=OP_arr_tr<tr_param_range(2,:) | OP_arr_tr>tr_param_range(1,:);
temp_result(:,[3 7])=0; % not examine CSF value
result(end+1,:)=temp_result;

result=sum(result(:));
result=~result;
end