
clc;

% load lookup table information
lkt_sim_set=load(fullfile('KB','sim_set.mat'));
lkt_sim_set=lkt_sim_set.sim_set;
lkt_layer_mus=load(fullfile('KB','layer_mus.mat'));
lkt_layer_mus=lkt_layer_mus.layer_mus;
lkt_mus_table=load(fullfile('KB','mus_table.txt'));

lkt_layer_mus={225,200,23,[25:25:275]}; %[100:25:200]
% lkt_layer_mus={[75:25:250],150,23,150};

in_place_arr=[];
% in_place_arr=zeros(length(lkt_layer_mus{1}),length(lkt_layer_mus{2}),length(lkt_layer_mus{3}),length(lkt_layer_mus{4}));
for L1=1:length(lkt_layer_mus{1})
    for L2=1:length(lkt_layer_mus{2})
        for L3=1:length(lkt_layer_mus{3})
            for L4=1:length(lkt_layer_mus{4})
                lkt_index=find(lkt_mus_table(:,1)==lkt_layer_mus{1}(L1) & lkt_mus_table(:,2)==lkt_layer_mus{2}(L2) & lkt_mus_table(:,3)==lkt_layer_mus{3}(L3) & lkt_mus_table(:,4)==lkt_layer_mus{4}(L4));
                in_place_arr(end+1)=lkt_index;
            end
        end
    end
end
    
mus_table=lkt_mus_table(in_place_arr,:);

dtof=[];
mua_param_arr=[0.45 0.3 0.042 0.4 0.2 0];
index=1;
for sim=in_place_arr
    load(fullfile('KB',['sim_' num2str(sim)],'PL_1.mat'));
    for s=1:5
        for g=1:10
            dtof(g,s,index)=1/each_photon_weight_arr(s)*sum(exp(-1*sum(double(SDS_detpt_arr{g,s}).*mua_param_arr,2)),1);
        end
    end
    index=index+1;
end


figure;
tiledlayout('flow');
for i=1:5
    nexttile;
    to_plot=squeeze(dtof(:,i,:));
    semilogy(1:1:10,to_plot);
%     ylim([min(to_plot(:)) max(to_plot(:))]);
end

relative_change=dtof(1:6,:,:)./dtof(1:6,:,6)-1;
figure;
tiledlayout('flow');
for i=1:5
    nexttile;
    to_plot=squeeze(relative_change(:,i,:));
    plot(1:1:6,to_plot);
    ylim([-1 0.5]);
    xlim([1 6]);
%     ylim([min(to_plot(:)) max(to_plot(:))]);
end