%{
Evaluate the influence of the error of OPs assumed to be known on interest
OPs

Ting-Yi Kuo
Last update: 2023/07/18
%}

clc;clear;close all;

global lambda Lbound Ubound net param_range;

%% param

Add_error_mode=0; % if =0, use result without error; =1 to use results with error

input_dir='test_fitting_2023-09-12-12-24-40'; % please move the fitting folders into this folder first.
subject_name_arr={'KB'};
num_anser_to_generate=25; % number of target spec (true answer)
num_error_to_generate=15; % number of adding noise to the same, the first one will have no error
times_to_fitting=20; % number of fitting using different init value

fitting_dir='fitting_SDS123456_345_gate1-5'; % the fitting folder
num_SDS_cw=6;
num_SDS_tr=5; % how many SDS are in the target spectrum
num_gate=10;
num_fitted_param=13; % the number of the fitted parameters

OP_total_result=zeros(2,4,length(fitting_dir));

if Add_error_mode==0
    to_process_fitting_index=1;
elseif Add_error_mode==1
    to_process_fitting_index=2:num_error_to_generate;
end


%% Arrange fitting OP result and target spec

OP_result=[];

for sbj=1:length(subject_name_arr)
    % make output folder and find fitting folder
    mkdir(fullfile(input_dir,'arrangement',subject_name_arr{sbj},fitting_dir));
    fitting_op_arrange=[];
    target_spec_arrange=[];

    for target_i=1:num_anser_to_generate
        fitting_op=[];
        for error_i=to_process_fitting_index

            temp_op_folder=fullfile(input_dir,subject_name_arr{sbj},['target_' num2str(target_i) '_' num2str(error_i)],fitting_dir);
            
            fitting_op_arrange_temp=load(fullfile(temp_op_folder,'fitting_result_sort.txt'));
            fitting_op(end+1,:)=fitting_op_arrange_temp(1,num_fitted_param+num_SDS_cw+num_SDS_tr*num_gate+3:end); %num_fitted_param+num_SDS+3:end

            temp_target_folder=fullfile(input_dir,subject_name_arr{sbj},['target_' num2str(target_i) '_' num2str(error_i)]);
            
            target_spec_arrange_temp=load(fullfile(temp_target_folder,'target_spec.txt'));
            if isempty(target_spec_arrange)
                target_spec_arrange(:,:,1,sbj)=target_spec_arrange_temp;
            else
                target_spec_arrange(:,:,end+1,sbj)=target_spec_arrange_temp;
            end
        end
        fitting_op_arrange(end+1,:)=mean(fitting_op,1);
    end
    OP_result(1:num_anser_to_generate,1:7,sbj)=fitting_op_arrange;
end

OP_result=mean(OP_result,3);

for i=1:3
    OP_result_arrange_temp(:,:,i)=OP_result(2+8*(i-1):1+8*i,:);
    OP_result_arrange(:,:,i)=[OP_result_arrange_temp(1:4,:,i); OP_result(1,:); OP_result_arrange_temp(5:8,:,i)];
end


%% Plot fitting OP result

xtitle_arr={'\mu_{s,scalp}','\mu_{s,skull}','\mu_{a,scalp}'};
figure('Position',[0 0 1920 480]);
ti=tiledlayout(1,3);
for i=1:3
    nexttile;
    for op=1:size(OP_result_arrange,2)-1
        p=plot(-4:1:4,100*OP_result_arrange(:,op,i),'-o');
        title(xtitle_arr{i});
        xticklabels({'-20%','-15%','-10%','-5%','0','5%','10%','15%','20%'});
        ylabel('error(%)')
        hold on;
    end
    p.Parent.YAxisLocation ="origin";
end
    
leg = legend('\mu_{a,scalp}','\mu_{s,scalp}','\mu_{a,skull}','\mu_{s,skull}','\mu_{a,GM}','\mu_{s,GM}','Orientation','horizontal');
leg.Layout.Tile = 'south';

mkdir('result');
if Add_error_mode
    print(fullfile(input_dir,'arrangement',subject_name_arr{sbj},fitting_dir,'OP_influence_error_Error.png'),'-dpng','-r200');
else
    print(fullfile(input_dir,'arrangement',subject_name_arr{sbj},fitting_dir,'OP_influence_error_noError.png'),'-dpng','-r200');
end


%% Calculate relative change

if Add_error_mode==0
    relative_change=zeros(num_gate,num_SDS_tr,num_anser_to_generate-1);

    for i=2:num_anser_to_generate
        relative_change(:,:,i-1)=target_spec_arrange(:,:,i)./target_spec_arrange(:,:,1)-1;
    end

    % compare change
    title_arr={'\mu_{s,scalp}','\mu_{s,skull}','\mu_{a,scalp}'};

    f=figure('Position',[0 0 1920 1080]);
    ti=tiledlayout(5,3);%(3,4)
    for s=1:num_SDS_tr
        for l=1:3
            nexttile;
            ind=8*(l-1)+1;
            for ch=1:8
                plot(1:1:num_gate,relative_change(:,s,ind));
                ind=ind+1;
                hold on;
            end
            xlabel('Time gate');
            ylabel('\DeltaR/R_{b}');
            xticks(1:10);
            xlim([1 10]);
            title(['SDS' num2str(s) ',' title_arr{l}]);
        end
    end
    leg = legend('-20%','-15%','-10%','-5%','5%','10%','15%','20%','Orientation','horizontal');
    leg.Layout.Tile = 'south';

    mkdir('results');
    print(fullfile('results','relative_change.png'),'-dpng','-r200');

end



