% OPs: 8 column of OPs, mua1 mus1 mua2 ......

function spec=fun_ANN_forward(OPs)
global net param_range;
% param
do_normalize=1; % if =1, do normalization to spec and param

%% main
OPs=OPs(:)';

if do_normalize
    param_scaling=param_range(1,:)-param_range(2,:);
    OPs=(OPs-param_range(2,[2 4 8]))./param_scaling([2 4 8]);
end
spec=predict(net,OPs');
spec=double(spec');
if do_normalize
    spec=power(10,-spec);
end

end
