%{
Find the range of K for each A, and construct a table for interpretation

Benjamin Kao
Last update: 2021/01/12
%}

clc;clear;close all;

%% param
subject_name_arr={'KB'}; % the name of the subjects
g=0.9; % g for L1,2,4
target_spec_wl=(700:900)';

model_folder='model_arrange'; % the folder of the models

%            A1     A2     A4      K1     K2     K4     
Lbound=     [5    0.1    0.1     1.0    1.0    1.0];
Ubound=     [35     150    150     1.5    1.5    5.0];

%% main
for sbj_i=1:length(subject_name_arr)
    % load the param range of the mus
    param_range=load(fullfile(model_folder,[subject_name_arr{sbj_i} '_cw_model.mat']), 'cw_param_range');
    param_range=param_range.cw_param_range(:,[5 6 8]); % only use mus1,2,4, (2,:) is min, (1,:) is max

    Lbound=[];
    Ubound=[];

    A_Krange_arr={};
    for L_i=1:size(param_range,2) % mus1,2,4
        Lbound(L_i*2-1)=param_range(2,L_i)*(1-g);
        Ubound(L_i*2-1)=param_range(1,L_i)*(1-g);
        A_interval=find_A_interval(Lbound(L_i*2-1),Ubound(L_i*2-1));
        K_range=[];
        for i=1:length(A_interval)
            K_range(i,:)=fitting_K_range(A_interval(i),target_spec_wl,g,param_range(2,L_i),param_range(1,L_i));
        end
        A_Krange_arr{L_i}=[A_interval K_range];

        Lbound(L_i*2)=min(K_range(:));
        Ubound(L_i*2)=max(K_range(:));
    end

    %% plot
    for L_i=1:size(param_range,2)
        nexttile();
        plot(A_Krange_arr{L_i}(:,1),A_Krange_arr{L_i}(:,2:3));
        legend({'K ub','K lb'},'Location','best');
        xlabel('A value');
    end

    saveas(gcf,fullfile(model_folder,['A_Krange_' subject_name_arr{sbj_i} '.png']));
    save(fullfile(model_folder,['A_Krange_arr_' subject_name_arr{sbj_i} '.mat']),'A_Krange_arr','Lbound','Ubound');
    close all;
end
disp('Done!');


%% functions
function output=fitting_K_range(A,wl,g,min_mus,max_mus)
% fitting the max K
K0=0;
output(1,1)=fmincon(@(x)-x,K0,[],[],[],[],[],[],@(x)AK_mus_constrain(x,A,g,wl,min_mus,max_mus));

% fitting the min K
output(1,2)=0;
end

function [c,ceq]=AK_mus_constrain(x,A,g,wl,min_mus,max_mus)
ceq=[];
mus=AK_to_mus(A,x,g,wl);
c=[max(mus)-max_mus ; min_mus-min(mus); -x];
end

function mus=AK_to_mus(A,K,g,wl)
mus = A * ( wl / 800 ) .^ -K ./ ( 1 - g );
end

% output the A array in the range, which is a geometric sequence (fix ratio)
function output=find_A_interval(min_rnage,max_range)
fix_ratio=1.05;
output=[];
now_value=min_rnage;
while now_value<=max_range
    output(end+1,1)=now_value;
    now_value=now_value*fix_ratio;
end
if max(output)~=max_range
    output(end+1)=max_range;
end
end