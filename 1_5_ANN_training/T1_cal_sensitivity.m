%{
Set the change rate for target layer to test the sensitivity of different SDS and gate

Ting-Yi Kuo
Last update: 2023/06/22
%}

close; clear all;
global net param_range;

baseline=[0.31 0.33 245];
changerate_to_exam=[-20 -10 10 20];
subject_name_arr={'KB'};
interpolation={'KB_test1_2023-07-10-18-40-48','BY_test1_2023-07-10-18-43-48'};
model_dir='model_arrange';
title_arr={'\mu_{a,skull}','\mu_{a,GM}','\mu_{s,GM}'};
ylim_arr={[-1 3],[-1 10],[-0.5 0.5]};
num_SDS=6;
num_gate=10;

%% make the mu table
mu_param_arr=[];
testing_index=1;

%% for baseline
temp_mu_arr=baseline;
temp_param=zeros(1,3);
temp_param(1,:)=temp_mu_arr;
mu_param_arr(testing_index,:)=temp_param;

testing_index=testing_index+1;

%% for changing mu
for l=1:size(baseline,2) % the layer to change
    for del_mus=1:length(changerate_to_exam)
        temp_mu_arr=baseline;
        temp_mu_arr(l)=temp_mu_arr(l)*(1+changerate_to_exam(del_mus)/100);
        temp_param=zeros(1,3);
        temp_param(1,:)=temp_mu_arr;
        mu_param_arr(testing_index,:)=temp_param;

        testing_index=testing_index+1;
    end
end

%% get dtof

dtof=zeros(num_gate,num_SDS,size(mu_param_arr,1));
temp_dtof=zeros(num_gate,num_SDS);

% ANN
% for sbj=1:length(subject_name_arr)
%     load(fullfile(model_dir,[subject_name_arr{sbj} '_model.mat'])); % net, param_range
%     for i=1:size(mu_param_arr,1)
%         temp_dtof_=fun_ANN_forward(mu_param_arr(i,:));
%         for s=1:num_SDS
%             temp_dtof(:,s)=temp_dtof_((s-1)*num_gate+1:(s-1)*num_gate+num_gate)';
%         end
%         dtof(:,:,i)=temp_dtof;
%     end
% end

% LUT
for sbj=1:length(subject_name_arr)
    load('/home/md703/Documents/Ty/1_4_MCX_lookup_table_dataGen/KB_test2_2023-06-26-23-37-43/all_param_arr.mat')
    for i=1:size(all_param_arr,1)
        temp_dtof_=all_param_arr(i,9:end);
        for s=1:num_SDS
            temp_dtof(:,s)=temp_dtof_((s-1)*num_gate+1:(s-1)*num_gate+num_gate)';
        end
        dtof(:,:,i)=temp_dtof;
    end
end

%% Calculate sensitivity

sens=zeros(num_gate,num_SDS,size(mu_param_arr,1)-1);

for i=2:9%size(mu_param_arr,1)
    sens(:,:,i-1)=dtof(:,:,i)./dtof(:,:,1)-1;
end

for i=10:13
    sens(:,:,i-1)=dtof(:,:,i)./dtof(:,:,14)-1;
end

% compare SDS
% f=figure('Position',[0 0 1800 600]);
% ti=tiledlayout(3,4);%(3,4)
% for l=1:size(mu_param_arr,1)-1
%     nexttile;
%     for s=1:4
%         plot(1:1:num_gate,sens(:,s,l)')
%         hold on;
%     end
% end
% legend('SDS1','SDS2','SDS3','SDS4','SDS5','SDS6');

% compare change
f=figure('Position',[0 0 1920 1080]);
ti=tiledlayout(6,3);%(3,4)
for s=1:num_SDS
    for l=1:size(baseline,2)
        nexttile;
        ind=4*(l-1)+1;
        for ch=1:4
            plot(1:1:num_gate,sens(:,s,ind));
            ind=ind+1;
            hold on;
        end
        xlabel('Time gate');
        ylabel('\DeltaR/R_{b}');
        title(['SDS' num2str(s) ',' title_arr{l}]);
        ylim(ylim_arr{l});
    end
end
leg = legend('-20%','-10%','10%','20%','Orientation','horizontal');
leg.Layout.Tile = 'south';
print('sensitivity_GM.png','-dpng','-r200');

figure;
semilogy(1:1:10,all_param_arr(10:13,39:48));
hold on;
semilogy(1:1:10,all_param_arr(1,39:48));
legend('-20%','-10%','10%','20%','baseline');
    
    


    