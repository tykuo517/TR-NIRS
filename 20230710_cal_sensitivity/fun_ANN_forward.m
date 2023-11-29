% OPs: 8 column of OPs, mua1 mus1 mua2 ......
% type: 0:continuous-wave, 1:time-resolved

function spec=fun_ANN_forward(OPs,type)
global cw_net cw_param_range tr_net tr_param_range;
% param
do_normalize=1; % if =1, do normalization to spec and param

%% main
% OPs=OPs(:,[1 3 5 7 2 4 6 8]);

% Choose type
if type==0
    net=cw_net;
    param_range=cw_param_range;
elseif type==1
    net=tr_net;
    param_range=tr_param_range;
    param_range=param_range(:,[1 2 4 5 6 8]);
%     OPs=OPs(:,[1 2 4 5 6 8]);
end

% Predict
if do_normalize
    param_scaling=param_range(1,:)-param_range(2,:);
    OPs=(OPs-param_range(2,:))./param_scaling;
end
spec=predict(net,OPs');
spec=double(spec');

if do_normalize
    spec=power(10,-spec);
end

end
