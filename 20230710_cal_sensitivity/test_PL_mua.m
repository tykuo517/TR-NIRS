

clear all;
load('/home/md703/Documents/Ty/20230704_testting_SDS/KB/sim_11/PL_1.mat')

target_mua=load('OP_sim_sen.txt');
target_mua=target_mua(:,1:4);
target_mua(:,5)=target_mua(:,4)*0.5;
target_mua(:,6)=0;

num_SDS=5;
num_gate=10;

dtof=zeros(size(target_mua,1),num_SDS*num_gate);
index=1;
for t=1:size(target_mua,1)
    index=1;
    for s=1:num_SDS
        for g=1:num_gate
            if size(SDS_detpt_arr{g,s},1)>0
                dtof(t,index)=1/each_photon_weight_arr(s)*sum(exp(-1*sum(double(SDS_detpt_arr{g,s}).*target_mua(t,:),2)),1);%*(true_r/sim_set.detector_r).^2;
            else
                dtof(t,index)=0;
            end
            index=index+1;
        end
    end
end

for i=2:size(dtof,1)
    dtof_change(i-1,:)=dtof(i,:); %./dtof(1,:)
end

for i=1:size(dtof_change,1)
    for s=1:num_SDS
        plot_pre(:,s,i)=dtof_change(i,1+num_gate*(s-1):num_gate*s);
    end
end

for s=1:num_SDS
    f=figure('Position',[0 0 1920 1080]);
    ti=tiledlayout("flow");
    for g=1:num_gate
        nexttile;
        plot(1:1:size(dtof_change,1),squeeze(plot_pre(g,s,:)));
%         xticks(1:9);
%         xticklabels({'-20', '-10', '0', '10', '20'});
%         xlabel('Time gate');
        ylabel('\DeltaR/R_{b}');
        title(['gate' num2str(g)]);
    end
%     xlabel('Time gate');
%     ylabel('\DeltaR/R_{b}');
%         xlim([1 10]);

end



figure;
for t=1:size(target_mua,1)
    semilogy(1:1:50,dtof(t,:));
    hold on;
    
end