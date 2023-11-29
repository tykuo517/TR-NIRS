%{
Return the constrain for the fitting function

Benjamin Kao
Last update: 2020/09/19
%}

function [c, ceq] = fun_mus_constrain(param_arr)
global scale_size A_Krange;

param_arr=param_arr.*scale_size;
ceq = [];

% exam the K for each A value
for L=1:3
    temp_Krange=interp1(A_Krange.A_Krange_arr{L}(:,1),A_Krange.A_Krange_arr{L}(:,2:3),param_arr(2*L-1),'pchip');
    if temp_Krange(1)<param_arr(2*L) || temp_Krange(2)>param_arr(2*L)
        c=1;
        return;
    end
end

[OP_arr,~]=fun_param_to_mu(param_arr,0);

% contraints matrix
if fun_in_OP_range(OP_arr)
    c=-1;
else
    c=1;
end

end