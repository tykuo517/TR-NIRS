%{
Auto fitting the spec using ANN
normalize the dimension of each AK param

Benjamin Kao
Last update: 2021/01/21
%}

function fitting_success=fun_auto_fittingMuaMus(param_init,output_dir)

%% param
global testing_index lambda output_folder final Lbound Ubound use_min_disk;
global target_spec orig_target_spec num_SDS_cw num_SDS_tr SDS_choosed_cw SDS_choosed_tr; % not used, but for save the fitting setting
output_folder=output_dir;
use_min_disk=1; % =1 to stop store too much files

global scale_size;
% %            A1     K1     A2      K2     A4     K4     hc_1    sto2_1 hc_2    sto2_2 hc_4     sto2_4  mel
% Lbound=     [0.1    0.1    0.1     1.0    1.0    1.0    1       0.3    1       0.3    50       0.3     0];
% Ubound=     [80     150    150     1.5    1.5    5.0    100     1      100     1      200      1       0.001];
scale_size=   [10     1      10      1      10     1      10      0.7    10      0.7    20       0.7     0.001];

%              A4     hc_2    sto2_2 hc_4     sto2_4
% scale_size=   [10     10      0.7    20       0.7   ];

%% normalize the scale of params
param_init=param_init./scale_size;
LLbound=Lbound./scale_size;
UUbound=Ubound./scale_size;

%% init
c=fun_mus_constrain(param_init);
if sum(find(c>0))>0
    fitting_success=0;
    return;
end

testing_index=1;
final=0;

assert(sum(param_init<=LLbound)==0,'initial value setting wrong!');
assert(sum(param_init>=UUbound)==0,'initial value setting wrong!');

save(fullfile(output_folder,'wl.txt'),'lambda','-ascii','-tabs');

%% fitting
global times;
times=0;
% options = optimoptions('fmincon','Algorithm','sqp','Display','iter','DiffMinChange',5*10^-5,'OptimalityTolerance',1e-7,'ConstraintTolerance',1e-9,'StepTolerance',1e-10,'MaxFunctionEvaluations',round(100*length(param_init)*1.5)); % increase the min step size for finding gradients
options = optimoptions('fmincon','Algorithm','sqp','Display','iter','DiffMinChange',5*10^-4,'OptimalityTolerance',1e-7,'ConstraintTolerance',1e-9,'StepTolerance',1e-10,'MaxFunctionEvaluations',round(100*length(param_init)*1.5)); % increase the min step size for finding gradients
param_final=fmincon(@fun_scale_param_error,param_init,[],[],[],[],LLbound,UUbound,@fun_mus_constrain,options);
final=1;
fun_forward_calError_chooseSDS(param_final.*scale_size);

if use_min_disk==0
    save(fullfile(output_folder,'setting_record.mat'));
end
fitting_success=1;

disp('Done!');
end

%% functions
function output=fun_scale_param_error(param_init)
global scale_size times;
output=fun_forward_calError_chooseSDS(param_init.*scale_size);
times=times+1;
end